# Golden Gate Assembly: Overhang Design Guide 

본 문서는 `overhang_automatation.py`에 구현된 **대규모(High-plex) 조립을 위한 오버행 설계 원칙**과 정밀 검증 알고리즘을 정리합니다.

---

## 1. 기본 설계 원칙 (MoClo Standard)

모든 후보 오버행은 다음 물리적/생물학적 제약을 기본적으로 통과해야 합니다.

- **연속 염기 제한**: 동일 염기 3개 연속 금지 (예: AAA, TTT x).
- **회문 구조 금지**: Self-ligation 방지를 위해 역상보 서열과 일치하는 서열 제외.
- **GC 함량 최적화**: 25% ~ 75% (1~3개) 유지로 결합 에너지 균형 확보.
- **제한효소/코돈 제거**: BsaI, BbsI 사이트 및 Start/Stop Codon 중복 금지.

---

## 2. 정밀 간섭 차단 알고리즘 (Balanced Fidelity)

30개 이상의 대규모 One-pot 조립을 위해 단순 해밍 거리를 넘어선 **고급 필터링**을 적용합니다.

### ① 상호 거리 및 역상보 검증 (HD $\ge 2$ & RC-HD $\ge 2$)
- 임의의 두 오버행 $S_1, S_2$는 최소 2개 이상의 염기가 달라야 합니다.
- 또한, 한 오버행과 다른 오버행의 **역상보 서열** 사이의 거리도 2 이상 유지하여 부품의 역방향 조립을 원천 차단합니다.

### ② 부분 일치 집중 차단 (3-Base Similarity Filter)
- 리가아제가 4bp 중 3bp만 맞아도 붙여버리는 실수를 방지하기 위해 다음 패턴을 추가로 검사합니다.
    - **N N _ N** 패턴 일치 금지 (앞 2글자와 마지막 1글자가 같음)
    - **N _ N N** 패턴 일치 금지 (첫 1글자와 뒤 2글자가 같음)

---

## 3. 대규모 조립(30개 이상) 전략

4bp 오버행 시스템의 수학적 한계($4^4=256$) 내에서 성공률을 유지하며 개수를 확보하는 전략입니다.

- **Greedy Search with Luck**: 무작위 샘플링과 검증을 수천 번 반복하여 최적의 조합을 찾아냅니다.
- **Middle-Match 완화**: 16개 이상의 세트가 필요할 경우 '가운데 2글자 중복 금지' 규칙을 '최소화' 수준으로 완화하여 목표 개수를 확보합니다.

---

## 4. Fidelity 계산 및 예측 모델

본 프로젝트의 `estimate_fidelity_strict` 함수는 NEBridge(NEB) 공식 도구의 데이터를 참고한 보수적 점수 체계입니다.

- **Score $\ge 90.0\%$**: 매우 안전. 10개 이상의 부품 조립에 권장.
- **Score $80 \sim 90\%$**: 실용적 범위. 실험 시 콜로니 확인 개수를 늘릴 것을 권장.
- **Score $< 80\%$**: 위험. 오버행 세트를 재설계하거나 계층적 조립(Hierarchical)을 고려해야 함.

---

## 5. 실험 추천 세트 (예시)

위 로직을 통해 추출된, Fidelity 90% 이상이 보장되는 5 세트입니다.
실제 NEBridge ligase fidelity에서 실행하여 추천로직을 검증하였습니다.


- 세트 1: GCCT, CTGC, ACCG, ACAC, TTGG (예상 성공률 100%)
    - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
    ![alt text](images/image.png)

    - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 99%
    ![alt text](images/image_1.png)

- 세트 2: AGTT, ACTG, GGAT, TCAT, GTTC (예상 성공률 100%)
   - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 99%
   - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
  
- 세트 3: TTGG, ACTG, GAGC, TTCA, GACT (예상 성공률 95%)
   - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
   - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
  
- 세트 4: AGTA, CTAT, TTCC, CCAC, CCTT (예상 성공률 95%)
  - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 98%
  - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
  
- 세트 5: GCGT, GTGG, GGTT, ACAG, TGGT (예상 성공률 100%)
  - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
  - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 100%
  

- 세트 6: CTGG, TTCG, TCAA, AGGC, CATT, GACG, AGCA, ATAC, GTTA, ATCT (예상 성공률 95%)
  - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 97%
  - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 98%

- 세트 7: GTAT, GTCA, TCGG, AGTG, CTTG, TCAA, GCTG, CGGT, TACA, AGGA (예상 성공률 95%)
  - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 99%
  - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 99%

- 세트 8: TGTA, ACGA, AGCG, TCTG, TTCG, GTGG, ATCC, CAAG, GACT, GCTC, CATT, TGCC, CTAT, TCAA, CACC, ACAT, CCTA, TTAC, CTGC, CTCA
  - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 93%
  - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 95%

- 세트 9 : ACAC,AACC,TACT,CTTA,CCGT,GCTC,CTAT,CAGG,ACTT,GGCA,CTGC,CAAC,GGAC,GAAT,CATT,CACA,ATCT,TCCG,CCAA,CTCG,GCCT,TTAC,TGGT,AGCG,CGTC,TCGC,TTCA,TCAT,TATC,TCTA
  - BsaI-HFv2, 37-16 cycling : Estimated ligation fidelity: 83%
  - BbsI-HFv2, 37-16 cycling : Estimated ligation fidelity: 87%
---

## 6. 대규모 조립 시 발생하는 문제점 (4bp 시스템의 한계)

BsaI, BbsI를 사용하는 표준 Golden Gate Assembly 시스템은 다음과 같은 물리적/수학적 한계에 부딪히게 됩니다.

- **수학적 설계 공간의 부족**: 4bp 오버행으로 만들 수 있는 모든 경우의 수는 $4^4 = 256$개뿐입니다. 여기서 안전 조건(GC 함량, 회문 금지 등)을 적용하면 실제 후보군은 100여 개로 급격히 줄어듭니다.
- **상호 간섭(Cross-talk) 확률 급증**: 조립하려는 부품 수가 많아질수록(예: 30개), 제한된 후보군 안에서 서로 '충분히 다른' 서열들을 고르는 것이 불가능해집니다.
- **Fidelity 하락의 딜레마**: 30개 부품을 한 번에 조립하려면 29개의 서로 다른 오버행이 필요합니다. 하지만 개수가 늘어날수록 오버행 간의 해밍 거리가 가까워져, 리가아제가 서열을 착각하고 엉뚱한 부품을 붙여버리는 **오조립(Mis-ligation)**이 발생할 확률이 매우 높아집니다.
- **결론**: 따라서 30개 이상의 부품을 One-pot으로 조립하면서 동시에 높은 성공률(Fidelity)을 확보하는 것은 기존 4bp 시스템만으로는 기술적 한계가 존재합니다.

## 7. 해결방안

30개 이상의 부품을 조립해야 할 경우, 4bp 시스템의 한계를 극복하기 위해 다음 전략들을 도입할 수 있습니다.

### ① 계층적 조립 (Hierarchical Assembly) - **가장 권장**
- **전략**: 30개를 한 번에 조립하는 대신, 5~6개씩 나누어 조립(Level 1)한 후, 완성된 중간 모듈들을 다시 합치는 방식(Level 2).
- **장점**: 각 단계에서 가장 Fidelity가 높은 검증된 오버행(Set 1~5)만 재사용할 수 있어 성공률이 비약적으로 상승합니다.

### ② PaqCI 효소 도입 (7bp Overhang)
- **전략**: BsaI(4bp) 대신 PaqCI(7bp) 제한효소를 사용합니다.
- **장점**: 가능한 조합의 수가 $4^7 = 16,384$개로 늘어납니다. 30개 부품을 조립하더라도 각 오버행 간의 해밍 거리를 3~4 이상으로 유지할 수 있어 One-pot 조립이 훨씬 안정적입니다.

### ③ 실험적 Fidelity Matrix 기반 정밀 선별
- **전략**: 단순 해밍 거리 계산이 아닌, NEB의 실제 리가아제 결합 데이터(Potapov Matrix)를 코드에 반영하여 선별합니다.
- **장점**: 이론적 거리보다 실제 실험에서 간섭이 적은 것으로 증명된 최적의 30개 조합을 소프트웨어적으로 찾아낼 수 있습니다.

### ④ 자동 조건 완화 알고리즘 (Progressive Relaxation)
- **전략**: 부품 수가 너무 많아(예: 25개 이상) 안전 조건을 모두 만족하는 오버행을 찾을 수 없을 때, 시스템이 자동으로 검증 강도를 낮추어 조립을 완수합니다.
- **단계**: Perfect (HD $\ge 2$ & RC $\ge 2$) $\rightarrow$ Good (HD $\ge 2$) $\rightarrow$ Minimal (HD $\ge 1$).
- **장점**: 설계 단계에서 오버행 부족으로 인해 워크플로우가 중단되는 것을 방지합니다.

### ⑤ 하이브리드 방식 (GGA + Gibson)
- **전략**: 15개는 Golden Gate로, 나머지 15개도 Golden Gate로 조립한 뒤, 이 두 덩어리를 Gibson Assembly로 합칩니다.
- **장점**: 큰 조각끼리 합칠 때는 4bp 오버행 제약이 없는 Gibson 방식을 써서 접합부 설계를 단순화합니다.


## 8. 향후 개발 방향
- module = 1 카셋으로 나눠서 조립한다면 많은 오버행이 필요하지 않을 것.
- Part - module 조립 후, multi module은 gibson assembly를 하는 방향으로

---
*Last Updated: 2026-02-25 (Optimization for 30-plex One-pot Assembly)*
