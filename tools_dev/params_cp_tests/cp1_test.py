import sys
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
from params import DesignParams, DEFAULT_PARAMS, PRIMER3_MAX_SIZE, THERMO_MAX_LEN

ok = True
def chk(label, cond, extra=""):
    global ok; ok = ok and cond
    print(f"  [{'OK ' if cond else 'FAIL'}] {label}" + (f" — {extra}" if extra else ""))

def raises(fn):
    try: fn(); return None
    except ValueError as e: return str(e)

print("1) 기본값 = 기존 하드코딩 값")
d = DEFAULT_PARAMS
chk("Tm 범위 55~65 (opt 60)", (d.tm_min, d.tm_max, d.tm_opt) == (55.0, 65.0, 60.0))
chk("길이 18~25 (opt 20)", (d.primer_len_min, d.primer_len_opt, d.primer_len_max) == (18,20,25))
chk("GC 20~80", (d.gc_min, d.gc_max) == (20.0,80.0))
chk("Tm 범위 폭 10", d.tm_max - d.tm_min == 10.0)
# 2026-11-05: Tm 축이 **NEB Tm Calculator 재현**으로 바뀌었다. NEB 버퍼는 Na+가 아니라
# K+/Tris 기반이라 na_mM=0이 정상이고, 기본값은 Q5 조건(0.5 µM = 500nM)이다.
# 근거는 basic/sequence_analysis.calculate_tm 독스트링 / params.NEB_TM_CONDITIONS.
chk("NEB Q5 버퍼 기본값 (Na 0 / K 50 / Tris 10 / Mg 1.5 / 500nM)",
    (d.na_mM, d.k_mM, d.tris_mM, d.mg_mM, d.primer_conc_nM) == (0.0, 50.0, 10.0, 1.5, 500.0))
chk("gibson overlap 30", d.gibson_overlap == 30)   # 2026-10-24 기본값 20 → 30 (연구실 지정)
chk("효소 BsaI", d.gg_enzyme == "BsaI")
# overhang/spacer는 필드가 아니라 **효소에서 파생**된다 (CP2)
chk("overhang/spacer는 효소 파생", (d.overhang_len, d.site_spacer_len) == (4,1),
    f"BsaI → overhang {d.overhang_len}, spacer {d.site_spacer_len}")
chk("BbsI는 spacer 2로 따라옴", d.with_(gg_enzyme="BbsI").site_spacer_len == 2)
# 2026-11-05: Tm 자를 primer3로 통일하면서 **축이 바뀐 값**이다(Biopython 55 ≈ primer3 59).
# 물리적 의미는 같다 — 근거는 basic/sequence_analysis.calculate_tm 독스트링.
chk("terminal Tm 59 (primer3 축)", d.terminal_tm_min == 59.0)
chk("poly-X 8", d.max_poly_x == 8)
chk("module 조립 linear", d.module_assembly_type == "linear")

print("\n2) 하드 리밋 — 에러 (실측 기반)")
e = raises(lambda: d.with_(primer_len_max=37).validate())
chk(f"primer_len_max 37 > primer3 내장 {PRIMER3_MAX_SIZE}", e is not None, (e or "")[:70])
#   2026-10-24(H1): `overlap + primer_len_min > THERMO_MAX_LEN`은 **하드 에러 → 경고**로
#   내려갔다. 그 전에는 "품질 계산이 아예 불가능"해서 막았지만, 지금은 >60bp를 ViennaRNA로
#   계산한다(경고 전용). 차단 게이트가 사라지는 건 사실이라 경고로 반드시 알린다.
#   ⚠️ overlap은 2026-10-25부터 **min–max 범위**다. min만 올리면 기본 max(45)보다 커져
#      "범위가 뒤집혔다"는 **다른** 에러로 걸린다 — 이 검사가 보려는 건 그게 아니므로
#      둘을 같이 올린다.
_hi = dict(gibson_overlap=50, gibson_overlap_max=50)
e = raises(lambda: d.with_(**_hi).validate())
chk("overlap 50은 더 이상 하드 에러가 아니다 (H1)", e is None, str(e)[:70])
_wl = d.with_(**_hi).validate() or []
chk(f"overlap 50 + min 18 > {THERMO_MAX_LEN} → 경고로 알린다",
    any("경고 전용" in w for w in _wl), "; ".join(_wl)[:70])
chk("기본값에는 이 경고가 없다(소음 방지)",
    not any("경고 전용" in w for w in (d.validate() or [])))
for label, kw in [("min>max", dict(primer_len_min=30)), ("opt 범위 밖", dict(primer_len_opt=40)),
                  # Na 0은 더 이상 오류가 아니다(NEB 버퍼는 K/Tris 기반) — 1가 염 **합산**이 0일 때만 막는다.
                  ("1가 염 전부 0", dict(na_mM=0, k_mM=0, tris_mM=0)),
                  ("GC min>max", dict(gc_min=90.0)),
                  ("미지원 효소", dict(gg_enzyme="EcoRI")),
                  ("assembly_type 오타", dict(module_assembly_type="ciruclar"))]:
    chk(label, raises(lambda k=kw: d.with_(**k).validate()) is not None)

print("\n3) 소프트 — 경고만, 차단 안 함")
# ⚠️ 2026-11-05(Y11) — 예전엔 78~88로 탐침했는데, 88℃가 이제 **입력 제한**(35~80)에
#   걸려 하드 에러가 된다. 소프트 경고(_SOFT_TM_RANGE 45~75)는 그보다 안쪽 띠에
#   살아 있으므로, 그 띠 안에서 탐침한다 — 입력 제한 안이면서 일반적 범위 밖인 창.
#   두 층이 각자 자기 일을 한다는 걸 같이 확인한다(소프트는 경고, 하드는 에러).
w = d.with_(tm_min=76.0, tm_max=80.0).validate()
chk("Tm 범위 76~80 → 경고 반환(에러 아님)", any("일반적 범위" in x for x in w), f"{len(w)}건")
chk("경고가 있어도 예외 없음", isinstance(w, list))
chk("입력 제한 밖(88℃)은 소프트가 아니라 하드 에러",
    "입력 허용 범위" in str(raises(lambda: d.with_(tm_min=78.0, tm_max=88.0).validate())))

print("\n4) canonical — 해시 재료")
chk("같은 값 → 같은 문자열", d.canonical() == DesignParams().canonical())
chk("다른 Tm → 다른 문자열", d.canonical() != d.with_(tm_max=68.0).canonical())
chk("float 자리수 고정", "tm_min=55.000000" in d.canonical())
# 필드명 정렬(사전순) — 새 필드가 늘어도 정렬 자체가 유지되는지로 검증
_keys=[x.split("=")[0] for x in d.canonical().split(";")]
chk("필드명 정렬", _keys == sorted(_keys), f"{_keys[0]} … {_keys[-1]}")

print("\n5) 직렬화 (galaxy 계약)")
chk("to_dict→from_dict 왕복 동일", DesignParams.from_dict(d.to_dict()) == d)
chk("일부 키만 있어도 복원", DesignParams.from_dict({"tm_max": 68.0}).tm_max == 68.0)
e = raises(lambda: DesignParams.from_dict({"tm_mx": 65.0}))
chk("모르는 키 → 조용히 무시하지 않고 에러", e is not None, (e or "")[:60])

print("\n6) 불변성")
try:
    d.tm_max = 70.0; chk("frozen 대입 차단", False)
except Exception:
    chk("frozen 대입 차단", True)
chk("with_()는 새 인스턴스", d.with_(tm_max=70.0) is not d and d.tm_max == 65.0)

print("\n=== CP1:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
