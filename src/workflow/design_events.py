"""설계 중 일어난 **우회·컷오프 사건**을 한곳에 모은다.

## 왜 필요한가

`build_log.txt`는 **백엔드가** 쓴다(`app/pipeline.py`의 `_Logger`). 엔진은 `print()`로만
남기고 stdout 리다이렉트가 없어서, 엔진이 실제로 한 우회 판단은 systemd journal로
흘러가고 **케이스 산출물에는 한 줄도 안 남았다**(2026-11-06 확인):

  - 고른 효소가 부품 내부에 있어서 다른 효소로 **전환**한 것
  - 접합부에 제한효소 사이트를 만드는 overhang을 **교체**한 것
  - 프라이머 설계가 안 돼서 anneal 시작을 **k nt 밀어** 재시도한 것
  - clean spacer 후보가 없어서 "가장 덜 접히는 것"으로 **내려앉은** 것

유저가 받는 건 "성공한 결과물"뿐이라, 그 결과가 **원래 설계대로인지 우회한 것인지**
구분할 방법이 없었다. 여기 모아서 `reports/design_events.csv`로 내보내고 백엔드가
build log에 이어 붙인다.

## 왜 모듈 전역인가

이벤트를 만드는 곳은 `assembly/`·`workflow/`·`basic/`에 흩어져 있고 호출 깊이가 5단을
넘는다(builder → prepare_parts → planner → primer_design). 기록 객체를 인자로 꿰면
그 경로의 함수 서명이 전부 바뀌므로 모듈 전역 하나로 둔다.

⚠️ **케이스마다 `reset()`을 불러야 한다.** 백엔드는 run 하나를 자식 프로세스 하나에서
돌리고 그 안에서 케이스를 **순차**로 돈다(`pipeline.py`의 `for ji, job in enumerate(jobs)`).
reset을 안 하면 case 2의 로그에 case 1의 사건이 섞인다.
"""
from dataclasses import dataclass, field
from threading import Lock
from typing import Any, Dict, List, Optional

#: 사건 종류 — CSV의 `kind` 열에 그대로 들어간다.
#: `detour`는 "설계를 바꿔서 통과시켰다", `limit`은 "한계에 걸려 품질을 낮췄다",
#: `blocked`는 "진행을 포기했다"다. 등급이 아니라 **성격** 구분이다.
KIND_DETOUR = "detour"
KIND_LIMIT = "limit"
KIND_BLOCKED = "blocked"
KIND_NOTE = "note"


@dataclass
class DesignEvent:
    """사건 하나.

    Attributes:
        stage: 어느 단계인지 (`part_prep` / `junction` / `spacer` / `primer` / `assembly`)
        kind: `KIND_*` 중 하나
        where: 어디서 (접합부 이름, 부품 이름, module scope 등)
        what: 무엇을 했는지 (사람이 읽는 한 줄, **영어**)
        why: 왜 그랬는지 — 걸린 기준·컷오프 값
        scope: 멀티모듈에서 어느 반응인지(`module 2/5`). 단일은 빈 문자열
    """

    stage: str
    kind: str
    where: str
    what: str
    why: str = ""
    scope: str = ""

    def as_row(self) -> Dict[str, Any]:
        return {"scope": self.scope, "stage": self.stage, "kind": self.kind,
                "where": self.where, "what": self.what, "why": self.why}

    def as_line(self) -> str:
        """build log에 붙일 한 줄."""
        mark = {KIND_DETOUR: "↪", KIND_LIMIT: "⚠", KIND_BLOCKED: "⛔"}.get(self.kind, "·")
        head = f"{mark} [{self.stage}] {self.where}: {self.what}"
        return f"{head} ({self.why})" if self.why else head


_lock = Lock()
_events: List[DesignEvent] = []
_scope = ""


def reset(scope: str = "") -> None:
    """기록을 비운다. **케이스를 시작할 때마다** 부른다."""
    global _scope
    with _lock:
        _events.clear()
        _scope = scope


def set_scope(scope: str) -> None:
    """이후 기록에 붙을 반응 범위(`module 2/5`). 멀티모듈 루프가 설정한다."""
    global _scope
    with _lock:
        _scope = scope


def record(stage: str, kind: str, where: str, what: str, why: str = "") -> None:
    """사건 하나를 남긴다.

    Args:
        stage: `part_prep` / `junction` / `spacer` / `primer` / `assembly`
        kind: `KIND_DETOUR` / `KIND_LIMIT` / `KIND_BLOCKED` / `KIND_NOTE`
        where: 접합부·부품 이름
        what: 한 줄 설명(영어 — 산출물 문구는 전부 영어다)
        why: 걸린 기준·컷오프
    """
    with _lock:
        _events.append(DesignEvent(stage=stage, kind=kind, where=where,
                                   what=what, why=why, scope=_scope))


def events() -> List[DesignEvent]:
    with _lock:
        return list(_events)


def rows() -> List[Dict[str, Any]]:
    return [e.as_row() for e in events()]


def lines() -> List[str]:
    return [e.as_line() for e in events()]


#: `reports/` 아래 기본 파일명. 백엔드가 케이스 출력으로 복사하고 build log에 이어 붙인다.
FILENAME = "design_events.csv"

FIELDS = ["scope", "stage", "kind", "where", "what", "why"]


def write_csv(out) -> int:
    """`reports/design_events.csv`로 쓴다.

    ⚠️ **사건이 0건이어도 쓴다.** "우회가 없었다"와 "기록이 안 돌았다"를 파일 유무로
    구분하면 `cross_checks.csv`에서 겪은 모호함을 그대로 되풀이한다(base.py 참고).
    0건이면 헤더 + `clean` 요약 행 하나를 쓴다.

    Returns:
        쓴 사건 수(요약 행은 안 센다)
    """
    import csv as _csv

    data = rows()
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = _csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        if not data:
            w.writerow({"scope": "", "stage": "(all)", "kind": "clean",
                        "where": "(none)",
                        "what": "No detour or cutoff was recorded during design.",
                        "why": "Everything was designed within the declared conditions."})
        else:
            w.writerows(data)
    return len(data)
