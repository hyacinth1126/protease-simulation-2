#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
방정식 계산 과정 설명 스크립트
prep_raw.csv에서 MM_calibration_equations.csv 방정식을 구하는 전체 과정을 단계별로 보여줍니다.
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

def exponential_association(t, F0, Fmax, k):
    """
    Exponential Association 모델
    F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
    """
    return F0 + (Fmax - F0) * (1 - np.exp(-k * t))


def read_raw_data(filename='prep_data/prep_raw.csv'):
    """Raw 데이터 읽기"""
    first_row_df = pd.read_csv(filename, header=None, nrows=1, sep='\t')
    concentration_row = first_row_df.iloc[0].values[1:]
    
    header_row_df = pd.read_csv(filename, header=None, skiprows=[0], nrows=1, sep='\t')
    header_names = header_row_df.iloc[0].values
    
    df = pd.read_csv(filename, header=None, skiprows=[0, 1], sep='\t')
    df.columns = header_names
    
    time_col = df.columns[0]
    times = pd.to_numeric(df[time_col].values, errors='coerce')
    
    data = {}
    i = 1
    conc_idx = 0
    
    while i < len(df.columns):
        if conc_idx < len(concentration_row):
            conc_value = float(concentration_row[conc_idx * 3])
        else:
            break
        
        conc_name = f"{conc_value} ug/mL"
        
        value_col_idx = i
        value_series = df.iloc[:, value_col_idx]
        valid_mask = (~pd.isna(value_series)) & (value_series.astype(str) != '') & (value_series.astype(str) != 'nan')
        valid_indices = np.where(valid_mask)[0]
        valid_times = pd.to_numeric(df.iloc[valid_indices, 0], errors='coerce').values
        valid_values = pd.to_numeric(df.iloc[valid_indices, value_col_idx], errors='coerce').values
        
        valid_mask2 = ~pd.isna(valid_values)
        valid_times = valid_times[valid_mask2]
        valid_values = valid_values[valid_mask2]
        
        if len(valid_times) > 0:
            data[conc_name] = {
                'time': valid_times,
                'value': valid_values,
                'concentration': conc_value
            }
        
        i += 3
        conc_idx += 1
    
    return data


def fit_equation_detailed(times, values, conc_name):
    """
    각 농도별로 방정식을 구하는 상세한 과정을 보여줍니다
    """
    print(f"\n{'='*70}")
    print(f"🔬 농도: {conc_name}")
    print(f"{'='*70}")
    
    # 1단계: 데이터 확인
    print(f"\n📊 1단계: Raw 데이터 확인")
    print(f"   시간 포인트 수: {len(times)}")
    print(f"   시간 범위: {times[0]:.1f} ~ {times[-1]:.1f} 분")
    print(f"   형광값 범위: {np.min(values):.2f} ~ {np.max(values):.2f}")
    print(f"   데이터:")
    for t, v in zip(times, values):
        print(f"      t={t:5.1f}분 → F(t)={v:10.2f}")
    
    # 2단계: 초기값 추정
    print(f"\n🔍 2단계: 초기 파라미터 추정")
    F0_init = values[0]
    Fmax_init = np.max(values)
    k_init = 0.1
    
    print(f"   F0 (초기값 추정) = {F0_init:.2f}")
    print(f"   Fmax (최대값 추정) = {Fmax_init:.2f}")
    print(f"   k (초기 추정) = {k_init:.4f}")
    
    # 3단계: 모델 정의
    print(f"\n📐 3단계: Exponential Association 모델")
    print(f"   모델: F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]")
    print(f"   파라미터: F0, Fmax, k")
    
    # 4단계: 최적화 피팅
    print(f"\n⚙️ 4단계: Non-linear Curve Fitting (최적화)")
    print(f"   방법: scipy.optimize.curve_fit (Levenberg-Marquardt 알고리즘)")
    print(f"   목표: 실측값과 모델 예측값의 차이(잔차 제곱합)를 최소화")
    
    try:
        popt, pcov = curve_fit(
            exponential_association, 
            times, values,
            p0=[F0_init, Fmax_init, k_init],
            bounds=([-1000, F0_init, 0.001], [Fmax_init, Fmax_init * 3, 10]),
            maxfev=5000
        )
        F0, Fmax, k = popt
        
        # 오차 추정
        perr = np.sqrt(np.diag(pcov))
        
        print(f"   ✅ 피팅 성공!")
        print(f"   최적 파라미터:")
        print(f"      F0 = {F0:.4f} ± {perr[0]:.4f}")
        print(f"      Fmax = {Fmax:.4f} ± {perr[1]:.4f}")
        print(f"      k (k_rate) = {k:.4f} ± {perr[2]:.4f}")
        
    except Exception as e:
        print(f"   ❌ 피팅 실패: {e}")
        F0, Fmax, k = F0_init, Fmax_init, k_init
        perr = [0, 0, 0]
    
    # 5단계: 피팅 품질 평가
    print(f"\n📈 5단계: 피팅 품질 평가")
    fit_values = exponential_association(times, F0, Fmax, k)
    residuals = values - fit_values
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((values - np.mean(values)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"   잔차 제곱합 (SS_res) = {ss_res:.2f}")
    print(f"   전체 제곱합 (SS_tot) = {ss_tot:.2f}")
    print(f"   R² (결정계수) = {r_squared:.4f}")
    print(f"   최대 잔차 = {np.max(np.abs(residuals)):.2f}")
    print(f"   평균 잔차 = {np.mean(np.abs(residuals)):.2f}")
    
    # 잔차 표시
    print(f"\n   잔차 (실측값 - 예측값):")
    for t, obs, fit, res in zip(times, values, fit_values, residuals):
        print(f"      t={t:5.1f}분: {obs:10.2f} - {fit:10.2f} = {res:8.2f}")
    
    # 6단계: MM 파라미터 계산
    print(f"\n🧮 6단계: Michaelis-Menten 파라미터 계산")
    Vmax = k * (Fmax - F0)
    Km = (Fmax - F0) / 2
    
    print(f"   Vmax = k * (Fmax - F0)")
    print(f"        = {k:.4f} * ({Fmax:.2f} - {F0:.2f})")
    print(f"        = {Vmax:.2f}")
    print(f"   Km = (Fmax - F0) / 2")
    print(f"      = ({Fmax:.2f} - {F0:.2f}) / 2")
    print(f"      = {Km:.4f}")
    
    # 7단계: 최종 방정식
    print(f"\n📝 7단계: 최종 방정식 생성")
    
    # 실제 방정식 (F0 포함)
    if abs(F0) < 0.01:
        equation = f"F(t) = 0.00 + ({Fmax:.2f}) * [1 - exp(-{k:.4f}*t)]"
    else:
        equation = f"F(t) = {F0:.2f} + ({Fmax - F0:.2f}) * [1 - exp(-{k:.4f}*t)]"
    
    print(f"   방정식: {equation}")
    
    # MM_calibration_equations.csv 형식
    print(f"\n💾 MM_calibration_equations.csv 형식:")
    print(f"   Concentration: {conc_name}")
    print(f"   F0: {F0:.1f}")
    print(f"   Fmax: {Fmax:.2f}")
    print(f"   Vmax: {Vmax:.1f}")
    print(f"   Km: {Km:.5f}")
    print(f"   k_rate: {k:.4f}")
    print(f"   Equation: {equation}")
    
    return {
        'F0': F0,
        'Fmax': Fmax,
        'k': k,
        'Vmax': Vmax,
        'Km': Km,
        'R_squared': r_squared,
        'equation': equation
    }


def main():
    """메인 함수"""
    print("=" * 70)
    print("📚 MM_calibration_equations.csv 방정식 계산 과정 설명")
    print("=" * 70)
    print("\n이 스크립트는 prep_raw.csv에서 각 농도별 시간 곡선을 피팅하여")
    print("MM_calibration_equations.csv의 방정식을 어떻게 구하는지 보여줍니다.\n")
    
    # Raw 데이터 읽기
    print("📂 Raw 데이터 읽는 중...")
    raw_data = read_raw_data('prep_data/prep_raw.csv')
    print(f"✅ {len(raw_data)}개 농도 조건 발견\n")
    
    # 각 농도별로 상세하게 보여주기
    results = {}
    for conc_name, data in sorted(raw_data.items(), key=lambda x: x[1]['concentration']):
        times = data['time']
        values = data['value']
        
        result = fit_equation_detailed(times, values, conc_name)
        results[conc_name] = result
    
    # 전체 요약
    print(f"\n{'='*70}")
    print("📋 전체 요약")
    print(f"{'='*70}")
    print(f"\n{'농도':<15} {'F0':>8} {'Fmax':>10} {'k_rate':>10} {'Vmax':>10} {'Km':>10} {'R²':>8}")
    print("-" * 70)
    
    for conc_name, result in results.items():
        print(f"{conc_name:<15} {result['F0']:>8.1f} {result['Fmax']:>10.2f} "
              f"{result['k']:>10.4f} {result['Vmax']:>10.1f} "
              f"{result['Km']:>10.5f} {result['R_squared']:>8.4f}")
    
    print(f"\n✨ 설명 완료!")
    print(f"\n핵심 요약:")
    print(f"  1. 각 농도별로 시간-형광값 데이터 추출")
    print(f"  2. Exponential Association 모델: F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]")
    print(f"  3. scipy.optimize.curve_fit으로 최적 파라미터 찾기")
    print(f"  4. Vmax = k * (Fmax - F0), Km = (Fmax - F0) / 2 계산")
    print(f"  5. 방정식 문자열 생성하여 CSV에 저장")


if __name__ == "__main__":
    main()

