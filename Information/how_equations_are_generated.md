# MM_calibration_equations.csv 방정식 생성 과정

## 📋 개요

`prep_raw.csv`에서 `MM_calibration_equations.csv`의 방정식이 생성되는 전체 과정을 설명합니다.

## 🔄 전체 프로세스

### 1단계: Raw Data 읽기 (`read_raw_data`)

```python
# prep_raw.csv 파일 구조:
# - 첫 번째 행: 농도 값들 (0.3125, 0.3125, 0.3125, 0.625, ...)
# - 두 번째 행: 헤더 (time_min, mean, SD, N, mean, SD, N, ...)
# - 세 번째 행부터: 실제 시간-형광값 데이터
```

각 농도별로 시간(time)과 형광값(value) 데이터를 추출합니다.

### 2단계: 각 농도별 시간 곡선 피팅 (`fit_time_course`)

각 농도 조건의 시간-형광값 데이터에 **Exponential Association 모델**을 피팅합니다.

#### 사용하는 수학 모델

```
F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
```

**파라미터 설명:**
- `F0`: 초기 형광값 (t=0일 때)
- `Fmax`: 최대 형광값 (포화 지점)
- `k`: 반응 속도 상수 (rate constant)

#### 피팅 과정

1. **초기값 추정:**
   ```python
   F0_init = values[0]  # 첫 번째 데이터 포인트
   Fmax_init = np.max(values)  # 최대값
   k_init = 0.1  # 초기 추정값
   ```

2. **Non-linear Curve Fitting:**
   ```python
   from scipy.optimize import curve_fit
   popt, pcov = curve_fit(
       exponential_association, 
       times, values,
       p0=[F0_init, Fmax_init, k_init],
       bounds=([-1000, F0_init, 0.001], [Fmax_init, Fmax_init * 3, 10])
   )
   F0, Fmax, k = popt
   ```
   
   `curve_fit` 함수가 최소제곱법을 사용하여 실험 데이터에 가장 잘 맞는 `F0`, `Fmax`, `k` 값을 찾습니다.

3. **Michaelis-Menten 파라미터 계산:**
   ```python
   Vmax = k * (Fmax - F0)  # 초기 반응 속도
   Km = (Fmax - F0) / 2    # 반속도 지점 근사
   ```

### 3단계: 방정식 문자열 생성

피팅된 파라미터로부터 방정식 문자열을 생성합니다:

```python
# F0가 거의 0인 경우 (실제로는 F0 = 0에 가까움)
if abs(F0) < 0.01:
    eq_str = f"F(t) = 0.00 + ({Fmax:.2f}) * [1 - exp(-{k:.4f}*t)]"
else:
    eq_str = f"F(t) = {F0:.2f} + ({Fmax:.2f}) * [1 - exp(-{k:.4f}*t)]"
```

**예시:**
```
F(t) = 0.00 + (9333.98) * [1 - exp(-2.0101*t)]
```

이 방정식은:
- 시간 `t`에 따른 형광값 `F(t)`를 계산합니다
- `F0 = 0.00`: 초기 형광값이 0
- `Fmax = 9333.98`: 최대 형광값
- `k = 2.0101`: 반응 속도 상수 (값이 클수록 빠르게 포화)

### 4단계: CSV 파일 저장

각 농도별로 다음 정보를 저장합니다:

| 컬럼 | 설명 |
|------|------|
| `Concentration` | 농도 이름 (예: "0.3125 ug/mL") |
| `F0` | 초기 형광값 |
| `Fmax` | 최대 형광값 |
| `Vmax` | 초기 반응 속도 (k * (Fmax - F0)) |
| `Km` | Michaelis 상수 근사값 |
| `k_rate` | 반응 속도 상수 (k) |
| `Equation` | 완성된 방정식 문자열 |

## 📊 실제 예시

### 입력 데이터 (prep_raw.csv)
```
시간(min)  0.3125 농도 형광값
0          0
1          7445.16
5          8383.01
10         8692.20
...
```

### 피팅 결과
```
F0 = 0.0
Fmax = 9333.98
k = 2.0101
Vmax = 9210.0 (계산값)
Km = 0.255 (계산값)
```

### 생성된 방정식
```
F(t) = 0.00 + (9333.98) * [1 - exp(-2.0101*t)]
```

## 🔬 수학적 배경

### Exponential Association 모델이란?

이 모델은 **일차 반응(first-order reaction)**을 따르는 과정을 설명합니다:
- 시간이 지날수록 형광값이 증가
- 일정 값(Fmax)에 점근적으로 접근 (포화)
- 반응 속도는 현재 상태에 비례

### 왜 이 모델을 사용하나?

1. **GraphPad Prism 표준**: 생화학 데이터 분석의 표준 모델
2. **생물학적 의미**: 효소 반응의 시간 경과를 잘 설명
3. **수학적 안정성**: 피팅이 안정적이고 수렴이 빠름

## 💻 코드 위치

주요 함수:
- `read_raw_data()`: `prep.py` 19-107줄
- `fit_time_course()`: `prep.py` 137-217줄
- `exponential_association()`: `prep.py` 110-115줄
- 방정식 생성: `prep.py` 467-489줄

## ✅ 검증

생성된 방정식을 검증하려면:
1. 방정식에 시간값(t)을 대입
2. 계산된 F(t) 값이 실제 실험 데이터와 비교
3. R² 값으로 피팅 품질 확인 (1에 가까울수록 좋음)


