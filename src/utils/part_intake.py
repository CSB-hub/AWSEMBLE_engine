"""유저 입력 → 바이오파츠 추출 (util_function_add).

플라스미드 합성/풀테스트와 **독립적인 유저 친화 입력 레이어**:
  - parts_from_sequences : 유저 서열 리스트 → Dseqrecord 부품 (DB 등록 불필요, 조립 바로)
  - parts_from_genbank   : genbank feature 기반 부품 추출 + DB 서열 매칭(재사용/미등록 표시)
  - find_part_by_sequence: 서열로 Bioparts 부품 검색

조립엔 Dseqrecord 부품을 그대로 쓴다(`build_plasmid(parts, vector, db)`의 parts).
※ io.py(DatabaseImport)를 수정하지 않기 위해 DB 서열매칭을 여기 독립 함수로 둔다
  (실행 중인 풀테스트와 파일 공유 회피 → 완전 독립). 추후 리팩 시
  DatabaseImport.find_part_by_sequence로 합칠 예정.
"""
from typing import List, Optional, Dict, Tuple

from pydna.dseqrecord import Dseqrecord

# DB 부품 테이블 (= part_type). backbone/vector/enzyme/nt/overhang/insert_dna 등
# 비-부품 테이블은 제외(부품 검색 대상 아님).
# 정의는 params.py 한 곳에 있다 — 여기서 따로 적지 않는다.
from params import MYSQL_PART_TABLES as _PART_TABLES_TUPLE
_PART_TABLES = list(_PART_TABLES_TUPLE)

# GenBank feature type(소문자) → DB part_type 매핑
FEATURE_TYPE_MAP = {
    "promoter": "promoter",
    "rbs": "rbs",
    "ribosome_binding_site": "rbs",
    "cds": "cds",
    "gene": "cds",
    "terminator": "terminator",
}


def find_part_by_sequence(db, sequence: str, part_type: Optional[str] = None) -> Optional[Tuple[str, str]]:
    """서열로 Bioparts 부품을 검색 → (part_id, part_type) 또는 None.

    Args:
        db: DatabaseImport (조회용; _execute만 사용 — io.py 미수정)
        sequence: 검색할 부품 서열
        part_type: 주어지면 그 테이블만, None이면 전 부품테이블(_PART_TABLES) 검색

    Returns:
        (part_id, part_type) 매칭 시 / 없으면 None

    Raises:
        ValueError: part_type가 부품테이블이 아닐 때(오타/비부품 타입 조기 검출) — fix6
    """
    seq = sequence.upper()
    # part_type 지정 시 그 테이블만. 알 수 없는 타입은 조용한 전수검색(의도 삼킴) 대신 에러로 드러낸다.
    if part_type is None:
        tables = _PART_TABLES
    elif part_type in _PART_TABLES:
        tables = [part_type]
    else:
        raise ValueError(f"알 수 없는 part_type: {part_type!r} (허용: {_PART_TABLES})")
    for t in tables:
        rows = db._execute("Bioparts", f"SELECT ID FROM `{t}` WHERE Sequence = %s", (seq,))
        if rows:
            return (rows[0][0], t)
    return None


def parts_from_sequences(
    sequences: List[str],
    types: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
) -> List[Dseqrecord]:
    """유저 서열 리스트 → Dseqrecord 부품 리스트 (조립 입력용, DB 무관).

    Args:
        sequences: 부품 서열 리스트 (조립 순서대로)
        types: 각 부품의 part_type (promoter/rbs/cds/terminator) — 선택, feature 기입용
        names: 각 부품 이름 — 선택 (없으면 user_part_N)

    Returns:
        Dseqrecord 리스트 (annotations['part_type'] 설정 시 module feature에 반영)
    """
    from basic.sequence_analysis import validate_dna_sequence

    n = len(sequences)
    # 병렬 리스트 길이 일치 검증 — 불일치 시 인덱싱 IndexError 대신 명확한 에러
    if types is not None and len(types) != n:
        raise ValueError(f"types 길이({len(types)})가 sequences 개수({n})와 다릅니다.")
    if names is not None and len(names) != n:
        raise ValueError(f"names 길이({len(names)})가 sequences 개수({n})와 다릅니다.")

    parts = []
    for i, seq in enumerate(sequences):
        s = str(seq).upper()
        # 입력 서열 조기 검증(ATGC-only; 빈 서열·N·IUPAC 거부) — 조립 단 E-1과 동일 게이트를
        # 입력 레이어에서 먼저 걸어 유저에게 명확한 에러를 준다.
        if not validate_dna_sequence(s):
            raise ValueError(
                f"부품 {i + 1} 서열이 유효하지 않습니다 "
                f"(ATGC만 허용, 빈 서열·N·IUPAC 불가): {seq!r}"
            )
        name = names[i] if names else f"user_part_{i + 1}"
        rec = Dseqrecord(s, name=name)
        if types:
            rec.annotations["part_type"] = types[i]
        parts.append(rec)
    return parts


def parts_from_genbank(gbk_path, db=None) -> List[Dict]:
    """genbank feature 기반 부품 추출 + DB 서열 매칭.

    각 feature(source 제외)를 part_type으로 매핑하고 서열을 추출한 뒤, DB에 같은 서열의
    부품이 있는지 확인한다(있으면 재사용, 없으면 신규로 표시).

    Args:
        gbk_path: genbank 파일 경로
        db: DatabaseImport (DB 매칭용; None이면 매칭 생략)

    Returns:
        리스트, 각 항목:
          {"part": Dseqrecord(annotations[part_type]),
           "part_type": str, "db_id": str|None,
           "status": "matched"(DB 재사용) | "new"(미등록) | "unknown_type"(매핑 안 됨),
           "label": str}
    """
    from utils.io import DNAFileReader

    # fix7: gbk는 여러 레코드(LOCUS)를 담을 수 있다. 현재는 첫 레코드만 처리하되,
    #       여러 개면 조용히 유실되지 않도록 경고한다(부품이 흩어져 있으면 분리 처리 필요).
    records = DNAFileReader.read_genbank(str(gbk_path))
    if len(records) > 1:
        print(f"⚠️ gbk에 레코드 {len(records)}개 발견 — 첫 번째만 처리합니다 "
              f"(나머지 {len(records) - 1}개 무시).")
    rec = records[0]

    # fix8: 부품 이름 유니크화 — 무명/동일 부품이 겹치면 build_plasmid registry(name 키)에서
    #       뒤엣것이 앞엣것을 덮어써 부품이 소실되므로, 겹치는 이름에 순번을 붙인다.
    _used: Dict[str, int] = {}

    def _uniq(base: str) -> str:
        n = _used.get(base, 0) + 1
        _used[base] = n
        return base if n == 1 else f"{base}_{n}"

    # 조립 부품의 코딩부는 CDS가 정답(실제 ORF)이고, gene은 같은 자리를 덮는 더 넓은 주석이다.
    # GenBank는 한 좌위에 gene+CDS를 함께 다는 관례가 있어 둘 다 추출하면 코딩부가 중복된다.
    # → CDS 구간을 미리 모아, gene은 그 자리에 CDS가 없을 때만 cds로 fallback(있으면 무시).
    cds_spans = [
        (int(f.location.start), int(f.location.end))
        for f in rec.features if f.type.lower() == "cds"
    ]

    out: List[Dict] = []
    for feat in rec.features:
        ftype = feat.type.lower()
        if ftype == "source":
            continue
        # gene: 같은 코딩영역을 CDS가 이미 커버하면 중복이므로 건너뛴다(CDS만 부품으로).
        if ftype == "gene":
            fs, fe = int(feat.location.start), int(feat.location.end)
            if any(fs < ce and cs < fe for (cs, ce) in cds_spans):
                continue
        ptype = FEATURE_TYPE_MAP.get(ftype)
        seq = str(feat.extract(rec).seq).upper()
        label = (feat.qualifiers.get("label") or feat.qualifiers.get("gene") or [""])[0]

        if ptype is None:
            out.append({
                "part": Dseqrecord(seq, name=_uniq(label or feat.type)),
                "part_type": feat.type, "db_id": None,
                "status": "unknown_type", "label": label,
            })
            continue

        match = find_part_by_sequence(db, seq, ptype) if db else None
        name = _uniq(match[0] if match else (label or f"new_{ptype}"))
        part = Dseqrecord(seq, name=name)
        part.annotations["part_type"] = ptype
        out.append({
            "part": part, "part_type": ptype,
            "db_id": match[0] if match else None,
            "status": "matched" if match else "new", "label": label,
        })
    return out


__all__ = ["parts_from_sequences", "parts_from_genbank", "find_part_by_sequence", "FEATURE_TYPE_MAP"]
