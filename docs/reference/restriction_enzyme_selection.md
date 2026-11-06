# 제한효소 선택 기준 및 출처 문헌

`src/basic/sequence_analysis.py`의 `find_restriction_sites()` 함수에 포함된
마스터 효소 목록(`main_enzymes`)의 선택 근거를 정리한 문서.

작성일: 2026-02-26

---

## 1. 효소 목록 전체 요약

| 효소 | 분류 | 선택 이유 |
|---|---|---|
| BsaI | Golden Gate Type IIS | 가장 표준적인 GGA 효소, 4bp overhang |
| BbsI | Golden Gate Type IIS | BsaI 대안, 동일한 4bp overhang |
| BsmBI | Golden Gate Type IIS | 레벨 2 조립에 표준으로 사용 |
| Esp3I | Golden Gate Type IIS | BsmBI 동의효소 |
| PaqCI | Golden Gate Type IIS | 7bp 인식 → unintended cut 확률 낮음, 대형 조립에 유리 |
| SapI | Golden Gate Type IIS | 7bp 인식, 3bp overhang, tRNA scaffold 클로닝 |
| EcoRI | 범용 | MCS 표준 |
| BamHI | 범용 | MCS 표준 |
| HindIII | 범용 | MCS 표준 |
| XhoI | 범용 | MCS 표준, Yeast HR 문헌 사용 예 |
| SalI | 범용 | XhoI compatible overhang |
| PstI | 범용 | MCS 표준 |
| XbaI | 범용 | BioBrick 표준 |
| NotI | 범용 | 8bp 인식, pET/pUC 계열 벡터 선형화 |
| SpeI | 범용 | BioBrick 표준, XbaI compatible overhang |
| NheI | 범용 | SpeI compatible overhang |
| AvrII | 범용 | NheI/SpeI compatible |
| KpnI | 범용 | pUC/pBR 계열 MCS 표준 |
| SmaI | 범용 | blunt-end cutter, Gibson 벡터 선형화 |
| ClaI | 범용 | pUC19 등 MCS 표준 |
| SphI | 범용 | Yeast HR 문헌 사용 예 |

---

## 2. Golden Gate Assembly용 Type IIS 효소

### 선택 기준
Type IIS 효소는 인식 서열 **외부**를 절단하여 인식 서열 흔적 없이 조립 가능.
4bp overhang을 생성하는 효소가 표준이며, NEB에서 권장하는 효소 우선순위를 따름.

### NEB 권장 우선순위
```
BsaI-HFv2 > PaqCI > BbsI-HF > SapI
```

### 효소별 세부 근거

**BsaI**
- Golden Gate Assembly에서 가장 널리 사용되는 표준 효소
- 6bp 인식 서열(`GGTCTC`), 4bp 5'-overhang 생성
- NEB에서 star activity를 줄인 BsaI-HF®v2 제공

**BbsI**
- BsaI를 사용할 수 없는 상황에서의 대안
- 6bp 인식 서열(`GAAGAC`), 4bp 5'-overhang 생성

**BsmBI / Esp3I**
- BsmBI와 Esp3I는 동의효소(isoschizomer)
- 계층적 조립(Hierarchical Assembly)에서 레벨 1과 레벨 2를 서로 다른 효소로 구분할 때 BsaI(레벨 1)와 쌍으로 사용

**PaqCI (AarI 동의효소)**
- 7bp 인식 서열(`CACCTGC`) → 인식 서열이 길어 part 내부에 우연히 존재할 확률이 낮음
- 50개 이상의 대형 fragment 조립에 유리
- AarI와 달리 star activity 없이 완전 절단 가능 (NEB 제공)
- 별도의 activator 시약 필요

**SapI (BspQI 동의효소)**
- 7bp 인식 서열, 3bp overhang 생성
- tRNA scaffold를 이용한 클로닝 및 특수 조립 프로토콜에 사용

---

## 3. Gibson Assembly / Yeast HR 클로닝용 범용 효소

### 선택 기준
Gibson Assembly와 Yeast Homologous Recombination(HR)은 서열 독립적(sequence-independent)이라
**어떤 효소든 벡터 선형화에 사용 가능**. 실제 문헌에서 자주 사용되는 효소들을 선정.

### 효소별 세부 근거

**NotI**
- 8bp 인식 서열(`GCGGCCGC`) → GC-rich, 인간/포유류 게놈에서 빈도 낮음
- pET, pUC 계열 등 광범위한 발현 벡터의 MCS에 포함

**SpeI / NheI / AvrII**
- BioBrick 표준(RFC 10)에서 XbaI/SpeI 쌍이 핵심 도구로 사용
- NheI, SpeI, AvrII, XbaI는 서로 compatible overhang(`CTAG`) 생성
  → 방향성 라이게이션 및 모듈식 조립에 유리

**KpnI / SmaI / ClaI**
- pUC19, pBR322 등 표준 클로닝 벡터의 MCS에 포함된 범용 효소

**SphI**
- Yeast HR 클로닝 문헌(Bhatt et al., 2014)에서 pET24a 벡터 선형화에
  `SphI + EcoRI + XhoI` 조합으로 직접 사용된 사례가 보고됨

---

## 4. BioPython 지원 여부 확인

아래 명령어로 전체 목록의 BioPython `Restriction` 모듈 지원 여부를 확인함.

```bash
conda run -n awsemble python -c "
from Bio.Restriction import Restriction
enzymes = ['BsaI','BbsI','BsmBI','Esp3I','PaqCI','SapI',
           'EcoRI','BamHI','HindIII','XhoI','SalI','PstI','XbaI',
           'NotI','SpeI','NheI','AvrII','KpnI','SmaI','ClaI','SphI']
[print(f'{e}: {hasattr(Restriction, e)}') for e in enzymes]
"
```

**결과: 21개 효소 전부 True** (2026-02-26 확인)

---

## 5. 출처 문헌

### Golden Gate Assembly

1. **NEB - Which restriction enzymes are used in Golden Gate Assembly?**
   https://www.neb.com/en-us/faqs/2017/07/17/which-restriction-enzymes-are-used-in-golden-gate-assembly

2. **NEB - Using PaqCI™ for Golden Gate Assembly**
   https://www.neb.com/en-us/nebinspired-blog/using-paqci-for-golden-gate-assembly

3. **Engler et al. (2008) - Golden Gate Shuffling: A One-Pot DNA Shuffling Method Based on Type IIs Restriction Enzymes**
   *PLOS ONE* 4(5): e5553
   https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0005553

4. **Marillonnet & Grützner (2020) - A User's Guide to Golden Gate Cloning Methods and Standards**
   *ACS Synthetic Biology* 12(1): 10-20
   https://pubs.acs.org/doi/10.1021/acssynbio.2c00355

### Gibson Assembly

5. **NEB Gibson Assembly® Master Mix Manual**
   https://www.neb.com/en/-/media/nebus/files/manuals/manuale2611.pdf

6. **Addgene - Plasmids 101: Gibson Assembly and Other Long-Homology Based Cloning Methods**
   https://blog.addgene.org/plasmids-101-gibson-assembly

### Yeast Homologous Recombination

7. **Bhatt et al. (2014) - A universal cloning method based on yeast homologous recombination that is simple, efficient, and versatile**
   *Microbial Cell Factories* 13: 120
   https://pmc.ncbi.nlm.nih.gov/articles/PMC4521215/

8. **Bhatt & Bhatt (2015) - Rapid and Efficient Plasmid Construction by Homologous Recombination in Yeast**
   *Cold Spring Harbor Protocols* 2015(9)
   https://cshprotocols.cshlp.org/content/2015/9/pdb.prot085100.abstract

### BioBrick 표준

9. **iGEM Registry - BioBrick RFC 10 Standard**
   https://parts.igem.org/Help:Standards/Assembly/RFC10
