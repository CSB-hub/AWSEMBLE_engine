"""
Assembly 공통 유틸리티 함수

DNA Assembly 전후 데이터 처리를 위한 공통 함수를 제공합니다.

핵심 기능:
    - remove_duplicate_features(): Assembly 결과물에서 중복 feature 제거
    - adjust_feature_positions(): 프라이머 추가로 인한 feature 위치 조정

참고:
    - KRIBB Line 130-170: remove_duplicate_features
    - KRIBB Line 497-504: adjust_feature_positions

작성일: 2025-11-06
"""

from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from Bio.Seq import Seq


def _circular_location(start: int, end: int, n: int, strand):
    """circular 서열(길이 n)에서 [start, end) 구간 location 생성 (origin wrap 처리).

    삽입 부품 feature가 plasmid 원점(origin)을 넘어가면 좌표를 0으로 되돌려(wrap)
    GenBank `join(start..n, 1..end')` 형식(CompoundLocation)으로 표현한다.
    (구버전은 경계를 넘는 feature를 그냥 버려 마지막 부품 feature가 소실됐다.)

    Args:
        start: 시작 좌표 (0-based, 0 이상)
        end: 끝 좌표 (exclusive, start 초과)
        n: circular 서열 길이
        strand: feature strand (+1/-1/None)

    Returns:
        FeatureLocation (경계 안) 또는 CompoundLocation (origin wrap)
    """
    if end <= n:
        # 경계 안: 단일 구간
        return FeatureLocation(start, end, strand=strand)
    if start >= n:
        # 전체가 경계 밖: 통째로 wrap
        return FeatureLocation(start - n, end - n, strand=strand)
    # origin을 가로지름: [start, n) + [0, end-n) 를 join
    parts = [
        FeatureLocation(start, n, strand=strand),
        FeatureLocation(0, end - n, strand=strand),
    ]
    if strand == -1:   # minus strand는 5'→3' 순서가 뒤집힘 (complement(join(...)))
        parts.reverse()
    return CompoundLocation(parts)


def annotate_inserted_region(
    construct: Dseqrecord,
    insert: Dseqrecord,
    target_start: int,
    target_end: int,
) -> Dseqrecord:
    """insert의 target 구간 feature를 construct의 삽입 위치에 매핑 기입 (P6).

    vector insertion(Gibson/GG) 후 조립 결과에는 vector feature만 남고 insert(부품)
    feature가 소실된다. construct에서 target 서열의 실제 위치를 찾아, target 구간 안에
    들어 있는 insert feature들을 그 위치로 좌표 변환(shift)해 복원한다.

    핵심 조립 로직과 무관한 **후처리**이며 서열을 바꾸지 않는다(SHA 불변).

    Args:
        construct: 조립 결과 (circular)
        insert: 원본 insert Dseqrecord (feature 보유)
        target_start: insert에서 증폭·삽입된 구간 시작 (0-based)
        target_end: 구간 끝 (exclusive)

    Returns:
        insert 부품 feature가 기입된 construct (제자리 수정 후 반환)
    """
    target_seq = str(insert.seq[target_start:target_end]).upper()
    if not target_seq:
        return construct

    cs = str(construct.seq).upper()
    n = len(cs)
    # circular 회전 대비 2배 문자열에서 위치 탐색 (정방향만 — vector insertion은 정방향)
    pos = (cs + cs).find(target_seq)
    if pos < 0:
        return construct
    pos %= n

    for f in insert.features:
        fs, fe = int(f.location.start), int(f.location.end)
        # target 구간 안에 완전히 포함되는 feature만 (부품 핵심 feature)
        if not (target_start <= fs < fe <= target_end):
            continue
        new_start = pos + (fs - target_start)
        new_end = pos + (fe - target_start)
        # circular 경계를 넘는 feature는 0으로 wrap → join() location으로 보존(구버전은 생략했음)
        construct.features.append(
            SeqFeature(
                _circular_location(new_start, new_end, n, f.location.strand),
                type=f.type,
                qualifiers=dict(f.qualifiers),
            )
        )
    return construct


# part_type(DB) → GenBank 표준 feature type. 미정의 type은 misc_feature로 둔다.
PART_TYPE_TO_FEATURE = {
    "promoter": "promoter",
    "rbs": "RBS",
    "cds": "CDS",
    "terminator": "terminator",
    # custom_part(유저가 융합해 둔 조합 파츠)는 GenBank 표준 type이 없다.
    # misc_feature로 두고 label(부품 id)로 식별한다 — 기본값과 같지만 **의도적**임을
    # 드러내려고 명시한다(빠뜨려서 misc_feature가 된 게 아니다).
    "custom_part": "misc_feature",
}


def annotate_module_parts(
    module: Dseqrecord,
    fragments: list,
    names: list,
    types: list = None,
    overhang_len: int = 4,
) -> Dseqrecord:
    """module에 각 part 구간 feature + 접합 overhang feature 기입 (P6, Golden Gate).

    Golden Gate 조립은 인접 fragment가 4bp overhang을 공유하며 이어지므로, fragment
    원본 통째로는(첫 부품 제외) module에 그대로 남지 않는다. 대신 앞 overhang을 뗀
    core(`frag[overhang_len:]`)는 보존되므로, 이를 exact match해 부품 위치를 잡는다.

    - part feature: 부품 본체 구간 [start, start+len]. feature type은 part_type을
      GenBank 표준(promoter/RBS/CDS/terminator)으로 매핑(미정의는 misc_feature),
      label=부품 id, note=part_type.
    - overhang feature: 인접 부품이 공유하는 4bp 접합 구간 (조립 구조 가시화)

    핵심 조립 로직과 무관한 **후처리**이며 서열을 바꾸지 않는다(SHA 불변).

    Args:
        module: 조립 결과 module
        fragments: 조립에 투입된 fragment(절단본) 리스트 (순서 = parts 순서)
        names: 각 fragment의 부품 id 리스트 (feature label)
        types: 각 부품의 part_type 리스트 (feature type 결정). None이면 전부 misc_feature.
        overhang_len: Golden Gate overhang 길이. **효소가 정하는 값**이므로 호출자가
            조립에 쓴 효소의 값(`GoldenGateModule.overhang_len`)을 넘겨야 한다.
            기본 4는 BsaI/BbsI/BsmBI 기준이며 SapI(3nt)에서는 틀린다.

    Returns:
        feature가 기입된 module (제자리 수정 후 반환)
    """
    ms = str(module.seq).upper()
    n = len(ms)
    # circular module(중간 벡터에 들어간 module)은 부품이 origin을 가로지를 수 있다.
    # 단순 find로는 그런 부품을 못 찾아 feature가 조용히 빠지므로, 2배 문자열에서 찾고
    # 좌표는 원래 길이로 되돌린다(위치 표현은 _circular_location이 join()으로 처리).
    is_circular = bool(getattr(module.seq, "circular", False))
    haystack = (ms + ms) if is_circular else ms

    if types is None:
        types = [""] * len(fragments)

    # 각 fragment core(앞 overhang 제외)로 exact match → part 시작 위치
    starts = []
    for f in fragments:
        fs = str(f.seq).upper()
        core = fs[overhang_len:] if len(fs) > overhang_len else fs
        p = haystack.find(core)
        starts.append((p - overhang_len) if p >= 0 else None)

    # part 본체 feature
    # extract_insert_fragment()가 만드는 fragment는 [left_oh][target][right_oh_rc] 구조라
    # len(f.seq)엔 양쪽 overhang이 다 포함돼있음 — 본체 feature는 그 core(target)만 가리켜야
    # 하므로 양 끝 overhang_len만큼 안쪽으로 당겨야 함. 단, 체인 맨 앞/뒤 fragment는
    # (S0002/S0005 같은 경계 spacer 포함) 그쪽 끝에 overhang이 아예 없는 blunt 구조라
    # (extract_insert_fragment의 site 유무 기반 front/back 트림 로직 참고) 그 바깥쪽은 안 당김.
    n_frags = len(fragments)
    for i, (f, name, ptype, st) in enumerate(zip(fragments, names, types, starts)):
        if st is None or st < 0:
            continue
        has_left_oh = i > 0
        has_right_oh = i < n_frags - 1
        start = st + (overhang_len if has_left_oh else 0)
        end = st + len(f.seq) - (overhang_len if has_right_oh else 0)
        if end <= start:
            continue
        if not is_circular and end > n:
            continue
        ftype = PART_TYPE_TO_FEATURE.get((ptype or "").lower(), "misc_feature")
        quals = {"label": [name]}
        if ptype:
            quals["note"] = [f"part_type={ptype}"]
        loc = (_circular_location(start % n, end, n, 1) if is_circular
               else FeatureLocation(start, end, strand=1))
        module.features.append(SeqFeature(loc, type=ftype, qualifiers=quals))

    # 인접 부품 공유 overhang feature (접합부 4bp)
    for i in range(len(fragments) - 1):
        st_i, st_j = starts[i], starts[i + 1]
        if st_i is None or st_j is None:
            continue
        end_i = st_i + len(fragments[i].seq)
        oh_start, oh_end = st_j, end_i   # 겹치는 overhang 구간
        if oh_end - oh_start == overhang_len:
            loc = (_circular_location(oh_start % n, oh_end, n, 1) if is_circular
                   else FeatureLocation(oh_start, oh_end, strand=1))
            module.features.append(
                SeqFeature(
                    loc,
                    type="misc_feature",
                    qualifiers={"label": ["overhang"], "note": [f"{names[i]}|{names[i+1]}"]},
                )
            )
    return module


def normalize_location(location):
    """길이 0인 조각을 떼어 낸 location을 돌려준다 (바뀔 게 없으면 원본 그대로).

    ## 왜 필요한가

    pydna가 환형 조립에서 `join{[5381:5381], [0:30]}`처럼 **염기 사이를 가리키는
    길이 0 조각**을 섞어 내보낸다. 서열상 의미는 `[0:30]` 하나와 같지만:

    - GenBank에 `join(5381^1,1..30)`으로 적히고, 숫자만 긁어 범위를 잡는 쪽(클라이언트
      벡터맵)이 `1..5381`로 읽어 **30bp feature가 지도를 통째로 덮는다**
    - 같은 feature가 `[0:30]`으로 한 번 더 들어와도 좌표가 달라 **중복으로 안 걸린다**

    둘 다 2026-10-25에 실제로 났다.

    Args:
        location: FeatureLocation 또는 CompoundLocation

    Returns:
        정규화된 location. 조각이 하나만 남으면 CompoundLocation을 풀어서 돌려준다.
    """
    parts = getattr(location, "parts", None)
    if not parts or len(parts) < 2:
        return location
    kept = [p for p in parts if int(p.end) > int(p.start)]
    if len(kept) == len(parts):
        return location
    if not kept:
        return location          # 전부 길이 0 — 건드리지 않는다(판단 보류)
    if len(kept) == 1:
        return kept[0]
    return CompoundLocation(kept)


def remove_duplicate_features(dseq: Dseqrecord) -> Dseqrecord:
    """Assembly 결과물에서 같은 위치의 중복 feature 제거
    # TODO: 현재 output gbk 파일에 feature를 기입하는 기능이 누락되어있음.
    # 위의 코드를 업데이트하면 이 코드는 자연스럽게 없어져도 되지않나? kribb 코드 참고했을때 만든 잔재인거 같음.

    Assembly 과정에서 vector와 insert 모두에 같은 feature가 있을 수 있습니다.
    예: vector와 insert 모두 "promoter" feature → 하나만 남김

    Args:
        dseq: Feature를 포함한 Dseqrecord

    Returns:
        중복이 제거된 Dseqrecord

    Notes:
        Feature identifier: (start, end, type, label)로 중복 판단

    Example:
        >>> # Assembly 결과물에 중복 feature가 있을 수 있음
        >>> result = gibson.assemble([vector, insert])[0]
        >>> result = remove_duplicate_features(result)
        >>> # 중복 제거됨

    참고:
        KRIBB modulation_kribb/code/assemble_pattern_modules.py Line 130-170
    """
    if not dseq.features:
        return dseq

    # Feature를 (start, end, type, label) 튜플로 변환하여 중복 체크
    seen = set()
    unique_features = []

    for feature in dseq.features:
        # 길이 0 조각을 먼저 떼어 낸다 — 안 그러면 같은 구간이 서로 다른 좌표로 보여
        # 중복 판정을 빠져나간다(§normalize_location).
        feature.location = normalize_location(feature.location)
        # Feature identifier 생성
        start = feature.location.start
        end = feature.location.end
        ftype = feature.type
        label = feature.qualifiers.get('label', [''])[0]

        identifier = (start, end, ftype, label)

        # 중복이 아니면 추가
        if identifier not in seen:
            seen.add(identifier)
            unique_features.append(feature)

    # 중복 제거된 features로 교체
    dseq.features = unique_features

    return dseq


def adjust_feature_positions(
    dseq: Dseqrecord,
    shift_length: int
) -> Dseqrecord:
    # TODO: remove_duplicate_features 함수와 같은 의문점.
    """프라이머 추가로 서열이 늘어난 경우 feature 위치 shift

    Gibson Assembly에서 프라이머로 5' 또는 3' 끝에 overlap을 추가하면
    모든 feature의 위치를 shift해야 합니다.

    예: 5' 끝에 20bp overlap 추가 → 모든 feature를 +20bp 이동

    Args:
        dseq: Feature를 포함한 Dseqrecord
        shift_length: 이동할 거리 (bp)
            - 양수: 오른쪽으로 이동 (5' 끝에 추가)
            - 음수: 왼쪽으로 이동 (잘라냄)

    Returns:
        위치가 조정된 Dseqrecord

    Example:
        >>> # 5' 끝에 20bp overlap 추가
        >>> dseq = adjust_feature_positions(dseq, shift_length=20)
        >>> # 모든 feature가 +20bp 이동됨

    참고:
        KRIBB modulation_kribb/code/assemble_pattern_modules.py Line 497-504
    """
    for feature in dseq.features:
        # 현재 위치
        start = feature.location.start
        end = feature.location.end
        strand = feature.location.strand

        # Shift 적용
        new_start = start + shift_length
        new_end = end + shift_length

        # 새 위치로 업데이트
        feature.location = FeatureLocation(new_start, new_end, strand=strand)

    return dseq


def annotate_module_parts_from_layout(module, records, log=None):
    """module에 부품 feature + 접합부 feature 기입 (homology 경로: Gibson / yeast HR).

    ## 왜 `annotate_module_parts`를 못 쓰나

    그쪽은 **Golden Gate 구조**를 가정한다 — fragment 하나 = 부품 하나, 인접 fragment가
    길이가 고정된 overhang(효소가 정하는 4nt)을 공유. homology 경로는 셋 다 다르다:

    1. 접합부 길이가 **자리마다 다르다**(품질을 통과하는 최소 길이를 골라서 40~60bp)
    2. arm보다 짧은 부품은 이웃 fragment에 흡수되므로 **fragment 하나가 부품 여러 개**를 담는다
    3. arm은 한쪽만 운반하므로 fragment 구조가 좌우 대칭이 아니다

    ## 그래서 탐색하지 않고 계산한다

    `JunctionPlan`이 배치를 이미 알고 있고 `FragmentSpec.to_record()`가 그걸 annotations에
    실어 보낸다(`part_spans` / `left_arm_len` / `body_len`). module 서열은
    `"".join(f.body for f in fragments)`이므로 **body 길이 누적합**이 곧 좌표다 —
    exact match로 되찾을 필요가 없고, 반복 서열이 있어도 틀리지 않는다.

    핵심 조립 로직과 무관한 **후처리**이며 서열을 바꾸지 않는다(SHA 불변).

    Args:
        module: 조립 결과 module (제자리 수정)
        records: 조립에 투입된 레코드 리스트 (순서 = fragment 순서).
            `FragmentSpec.to_record()`가 만든 것이어야 배치 정보가 실려 있다.
        log: 한 줄짜리 메시지를 받을 함수(선택). 없으면 `print`.

    Returns:
        feature가 기입된 module. 배치 정보가 없으면 그대로 돌려준다.
    """
    say = log or (lambda m: print(m))
    n = len(module.seq)
    is_circular = bool(getattr(module.seq, "circular", False))

    layouts = [r.annotations.get("part_spans") for r in records]
    if not any(layouts):
        say("  ⚠️  part feature 기입 생략: 배치 정보(part_spans)가 없습니다")
        return module

    # body 누적합이 module 좌표. 총합이 module 길이와 다르면 전제가 깨진 것이므로
    # **조용히 틀린 좌표를 쓰지 않고** 기입을 포기한다.
    total = sum(int(r.annotations.get("body_len") or 0) for r in records)
    if total != n:
        say(f"  ⚠️  part feature 기입 생략: body 합계 {total}bp != module {n}bp")
        return module

    offset = 0
    bounds = []
    for rec, spans in zip(records, layouts):
        for pid, ptype, a, b in (spans or []):
            start, end = offset + a, offset + b
            if end <= start or end > n:
                continue
            ftype = PART_TYPE_TO_FEATURE.get((ptype or "").lower(), "misc_feature")
            quals = {"label": [pid]}
            if ptype:
                quals["note"] = [f"part_type={ptype}"]
            loc = (_circular_location(start % n, end, n, 1) if is_circular
                   else FeatureLocation(start, end, strand=1))
            module.features.append(SeqFeature(loc, type=ftype, qualifiers=quals))
        offset += int(rec.annotations.get("body_len") or 0)
        bounds.append(offset)

    # 접합부 feature — 인접 fragment가 **공유하는** 구간. GG의 overhang에 대응하지만
    # 길이가 자리마다 다르고, **어느 쪽이 운반하느냐에 따라 경계의 좌/우가 갈린다**:
    #
    #   · 오른쪽이 운반(`left_arm`) → 그 fragment가 앞 fragment의 body **끝**까지 되짚는다
    #                                  → 창은 경계의 **왼쪽**
    #   · 왼쪽이 운반(`right_arm`)  → 그 fragment가 뒤 fragment의 body **앞**까지 넘어간다
    #                                  → 창은 경계의 **오른쪽**
    #
    # 한쪽으로 고정하면 절반이 엉뚱한 자리를 가리킨다(실측으로 잡은 버그).
    for i in range(len(records) - 1):
        left_carry = int(records[i].annotations.get("right_arm_len") or 0)
        right_carry = int(records[i + 1].annotations.get("left_arm_len") or 0)
        edge = bounds[i]
        if right_carry > 0:
            arm, start, end = right_carry, edge - right_carry, edge
        elif left_carry > 0:
            arm, start, end = left_carry, edge, edge + left_carry
        else:
            continue
        if start < 0 or end > n:
            continue
        loc = (_circular_location(start % n, end, n, 1) if is_circular
               else FeatureLocation(start, end, strand=1))
        module.features.append(SeqFeature(
            loc, type="misc_feature",
            qualifiers={"label": ["homology"],
                        "note": [f"{records[i].name}|{records[i+1].name} {arm}bp"]}))
    return module


# 모듈 메타데이터
__all__ = [
    'annotate_inserted_region',
    'normalize_location',
    'annotate_module_parts',
    'annotate_module_parts_from_layout',
    'remove_duplicate_features',
    'adjust_feature_positions',
]
