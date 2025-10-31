# 방정식 계산 과정 설명

## 개요

`prep_raw.csv`에서 각 농도별 시간 경과 곡선(time course)을 **Exponential Association 모델**로 피팅하여 방정식을 구합니다.

## 단계별 과정

### 1단계: Raw 데이터 읽기

`prep_raw.csv`에서 각 농도별로 시간(time_min)과 형광값(mean) 데이터를 추출합니다.

예시 데이터:
- 농도: 0.3125, 0.625, 1.25, 2.5, 5 ug/mL
- 각 농도마다 시간별 형광값 측정: 0분, 1분, 5분, 10분, 15분, 20분, 25분, 30분

### 2단계: 각 농도별 시간 곡선 피팅

각 농도에 대해 **Exponential Association 모델**을 피팅합니다:

```
F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
```

여기서:
- **F0**: 초기 형광값 (t=0일 때)
- **Fmax**: 최대 형광값 (t→∞일 때)
- **k**: 반응 속도 상수 (k_rate)
- **t**: 시간 (분)

### 3단계: 파라미터 최적화 (Non-linear Curve Fitting)

`scipy.optimize.curve_fit` 함수를 사용하여 실제 데이터에 가장 잘 맞는 F0, Fmax, k 값을 찾습니다.

```python
from scipy.optimize import curve_fit

# Exponential Association 모델 함수
def exponential_association(t, F0, Fmax, k):
    return F0 + (Fmax - F0) * (1 - np.exp(-k * t))

# 초기값 추정
F0_init = values[0]  # 첫 번째 데이터 포인트
Fmax_init = np.max(values)  # 최대값
k_init = 0.1  # 초기 추정값

# 최적화 피팅
popt, pcov = curve_fit(
    exponential_association, 
    times, values,
    p0=[F0_init, Fmax_init, k_init],
    bounds=([-1000, F0_init, 0.001], [Fmax_init, Fmax_init * 3, 10])
)

F0, Fmax, k = popt
```

### 4단계: MM 파라미터 계산

피팅된 값으로부터 Michaelis-Menten 파라미터를 계산합니다:

```python
# Vmax: 초기 반응 속도
Vmax = k * (Fmax - F0)

# Km: 반속도 지점 근사값
Km = (Fmax - F0) / 2
```

### 5단계: 방정식 생성

최종 방정식 문자열을 생성합니다:

```python
# 실제로는 F0가 거의 0이면:
equation = f"F(t) = {F0:.2f} + ({Fmax:.2f}) * [1 - exp(-{k:.4f}*t)]"

# MM_calibration_equations.csv 형식:
# F(t) = 0.00 + (9333.98) * [1 - exp(-2.0101*t)]
```

## 예시: 0.3125 ug/mL 농도

Raw data에서:
- 시간: 0, 1, 5, 10, 15, 20, 25, 30 (분)
- 형광값: 0, 7445.16, 8383.01, 8692.2, 9333.98, 9165.91, 9165.72, 9301.82

피팅 결과:
- F0 = 0.0 (초기값)
- Fmax = 9333.98 (최대 형광값)
- k_rate = 2.0101 (반응 속도 상수)

최종 방정식:
```
F(t) = 0.00 + (9333.98) * [1 - exp(-2.0101*t)]
```

## 전체 과정 요약

```
prep_raw.csv (시간별 형광값)
    ↓
각 농도별로 분리
    ↓
Exponential Association 모델 피팅
F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
    ↓
최적 파라미터 찾기 (curve_fit)
    ↓
Vmax, Km 계산
    ↓
방정식 문자열 생성
    ↓
MM_calibration_equations.csv 저장
```

## 중요한 포인트

1. **Non-linear curve fitting**: 실제 데이터를 수학적 모델에 맞추는 최적화 과정
2. **Exponential Association 모델**: 시간에 따라 지수적으로 증가하다가 Fmax에 수렴하는 형태
3. **Least squares method**: 실측값과 모델 예측값의 차이(잔차)를 최소화하는 파라미터를 찾음
4. **R² 값**: 피팅의 품질을 나타내는 지표 (1에 가까울수록 좋음)

