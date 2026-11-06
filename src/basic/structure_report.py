"""프라이머·결합부위 2차구조 그림 리포트 (HTML + 인라인 SVG).

## 왜 그림인가

실험자는 ΔG 숫자가 아니라 **구조를 보고** 판단한다 — "결합부가 뚫려 있으면 괜찮다"는
식이다(연구실 실제 작업 방식). 지금까지는 그 판단을 하려고 Benchling 같은 외부 도구를
따로 열어야 했다.

ΔG 경고를 화면에서 뺀 자리(2026-10-20)를 이 리포트가 메운다. 숫자는 표에 남고,
판단 근거인 그림은 여기 있다.

## 두 분자를 모두 그린다

  · **template 결합부위** — 프라이머가 앉을 자리가 접혀 막혔나.
    파츠 서열은 못 바꾸므로 조건(온도·첨가제)으로만 대응한다.
  · **프라이머 자체** — 프라이머가 자기 자신으로 hairpin을 만드나.
    이건 다시 설계하면 풀린다.

같은 "2차구조"라도 대응이 달라서 나눠 보여준다.
"""
import html
import os
import tempfile
from typing import Dict, List, Optional


#: viewBox를 내용 바깥으로 넓힐 여백(사용자 좌표). 글자는 **앵커를 중심으로** 그려지므로
#  앵커가 박스 안에 있어도 글리프 절반이 밖으로 나간다. 염기 하나가 차지하는 폭보다
#  넉넉히 잡는다.
_VIEWBOX_PAD = 40.0


def _fit_viewbox(svg: str) -> str:
    """ViennaRNA SVG에 `viewBox`를 넣고 고정 크기를 떼어낸다.

    ## 왜 필요한가

    `RNA.svg_rna_plot`이 내보내는 SVG는 `width="452" height="452"`만 있고
    **viewBox가 없다.** 그래서:

    1. CSS로 `max-width:100%`를 줘도 **좌표계가 같이 줄지 않아** 폭이 좁아지면
       그림이 그대로 잘린다(팝업처럼 좁은 칸에서 바로 드러난다).
    2. 염기 글자는 앵커를 중심으로 그려지는데 앵커가 박스 가장자리(실측 x≈446/452)까지
       가므로, 452 안에 "들어 있어도" 글리프 절반이 밖으로 나가 잘린다.

    내용의 실제 경계를 구해 여백을 더한 viewBox로 바꾸면 둘 다 해결된다.
    좌표는 `<g transform="scale(sx,sy) translate(tx,ty)">` 안에 있으므로 변환을 적용해
    계산한다 — 안 그러면 엉뚱한 범위가 나온다.

    Args:
        svg: ViennaRNA가 만든 SVG 문자열

    Returns:
        viewBox가 들어간 SVG. 파싱에 실패하면 원본을 그대로 돌려준다(그림이 없는 것보다 낫다).
    """
    import re

    m = re.search(r"<svg[^>]*?>", svg)
    if not m:
        return svg
    head = m.group()
    dim = {k: float(v) for k, v in re.findall(r'\b(width|height)="([\d.]+)"', head)}
    w, h = dim.get("width", 0.0), dim.get("height", 0.0)
    if w <= 0 or h <= 0:
        return svg

    g = re.search(r'transform="scale\(([\d.]+),([\d.]+)\)\s*translate\(([-\d.]+),([-\d.]+)\)"', svg)
    sx, sy, tx, ty = (float(g.group(i)) for i in (1, 2, 3, 4)) if g else (1.0, 1.0, 0.0, 0.0)

    xs, ys = [], []
    for pat, into in ((r'\bx="([-\d.]+)"', xs), (r'\bx1="([-\d.]+)"', xs), (r'\bx2="([-\d.]+)"', xs),
                      (r'\by="([-\d.]+)"', ys), (r'\by1="([-\d.]+)"', ys), (r'\by2="([-\d.]+)"', ys)):
        into.extend(float(v) for v in re.findall(pat, svg))
    # 배경 <rect>의 0 / w / h는 변환 밖 좌표라 섞이면 범위가 틀어진다.
    xs = [v for v in xs if v not in (0.0, w)]
    ys = [v for v in ys if v not in (0.0, h)]
    if not xs or not ys:
        return svg

    rx = [(v + tx) * sx for v in xs]
    ry = [(v + ty) * sy for v in ys]
    # 선언된 박스와 실제 내용의 **합집합**을 쓴다 — 내용이 박스 밖으로 나가는 경우도 있다.
    x0 = min(0.0, min(rx)) - _VIEWBOX_PAD
    y0 = min(0.0, min(ry)) - _VIEWBOX_PAD
    x1 = max(w, max(rx)) + _VIEWBOX_PAD
    y1 = max(h, max(ry)) + _VIEWBOX_PAD

    new_head = re.sub(r'\s*\b(width|height)="[\d.]+"', "", head)
    new_head = new_head.replace(
        "<svg",
        f'<svg viewBox="{x0:.0f} {y0:.0f} {x1 - x0:.0f} {y1 - y0:.0f}" '
        f'preserveAspectRatio="xMidYMid meet"',
        1,
    )
    return svg.replace(head, new_head, 1)


def _svg(seq: str, structure: str) -> str:
    """ViennaRNA로 구조 SVG를 그려 문자열로 돌려준다. 실패하면 빈 문자열."""
    try:
        import RNA
    except ImportError:
        return ""
    tmp = tempfile.NamedTemporaryFile(suffix=".svg", delete=False)
    tmp.close()
    try:
        RNA.svg_rna_plot(seq, structure, tmp.name)
        with open(tmp.name, encoding="utf-8") as f:
            svg = f.read()
        # 인라인 삽입을 위해 XML 선언 제거
        svg = svg.split("?>", 1)[-1] if "?>" in svg else svg
        return _fit_viewbox(svg)
    except Exception:
        return ""
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass


def fold(seq: str):
    """(구조 dot-bracket, MFE) — ViennaRNA 없으면 (None, None)."""
    try:
        import RNA
    except ImportError:
        return None, None
    try:
        ss, mfe = RNA.fold_compound(seq.upper()).mfe()
        return ss, mfe
    except Exception:
        return None, None


def _panel(title: str, seq: str, note: str = "") -> str:
    ss, mfe = fold(seq)
    if ss is None:
        return f"<div class='p'><h4>{html.escape(title)}</h4><p>ViennaRNA not available</p></div>"
    svg = _svg(seq, ss)
    return (
        f"<div class='p'><h4>{html.escape(title)}</h4>"
        f"<p class='m'>MFE {mfe:.2f} kcal/mol · {len(seq)} nt{' · ' + html.escape(note) if note else ''}</p>"
        f"<pre>{html.escape(seq)}\n{html.escape(ss)}</pre>"
        f"<div class='svg'>{svg}</div></div>"
    )


def build_html(fragments: List[Dict]) -> str:
    """fragment별 구조 리포트 HTML.

    fragments: {name, fwd_primer, rev_primer, fwd_site, rev_site,
                fwd_unpaired, rev_unpaired, notes} 목록
    """
    css = """body{font-family:system-ui,-apple-system,sans-serif;margin:24px;color:#243b53}
h2{border-bottom:2px solid #2c6ba0;padding-bottom:6px;margin-top:36px}
.row{display:flex;flex-wrap:wrap;gap:18px}
.p{flex:1 1 320px;border:1px solid #dfe4e6;border-radius:8px;padding:12px;background:#fff}
.p h4{margin:0 0 4px;font-size:.92rem}
.m{margin:0 0 8px;font-size:.78rem;color:#6b7a88}
pre{font-size:.68rem;overflow-x:auto;background:#f7f9fa;padding:8px;border-radius:4px;margin:0 0 8px}
.svg svg{width:100%;height:auto;max-height:420px;display:block;background:#fff}
.note{background:#fdf6e3;border-left:3px solid #d9a441;padding:8px 12px;font-size:.8rem;margin:8px 0}
.lead{color:#5b6b7a;font-size:.86rem;max-width:760px}"""
    out = [
        "<!doctype html><html lang='ko'><head><meta charset='utf-8'>",
        "<title>Secondary structure report</title>",
        f"<style>{css}</style></head><body>",
        "<h1>Secondary structure report</h1>",
        "<p class='lead'>The <b>template binding site</b> a primer sits on and the "
        "<b>primer itself</b> are folded separately. A blocked binding site cannot be "
        "fixed by changing the part sequence - handle it with conditions (temperature, "
        "additives). Structure in the primer itself is solved by redesign.</p>",
    ]
    for f in fragments:
        out.append(f"<h2>{html.escape(str(f.get('name', '')))}</h2>")
        for n in f.get("notes") or []:
            out.append(f"<div class='note'>{html.escape(str(n))}</div>")
        out.append("<div class='row'>")
        for label, seq, up in (
            ("Fwd binding site (template)", f.get("fwd_site"), f.get("fwd_unpaired")),
            ("Rev binding site (template)", f.get("rev_site"), f.get("rev_unpaired")),
            ("Fwd primer", f.get("fwd_primer"), None),
            ("Rev primer", f.get("rev_primer"), None),
        ):
            if not seq:
                continue
            note = f"unpaired probability {up:.2f}" if isinstance(up, (int, float)) else ""
            out.append(_panel(label, str(seq), note))
        out.append("</div>")
    out.append("</body></html>")
    return "\n".join(out)
