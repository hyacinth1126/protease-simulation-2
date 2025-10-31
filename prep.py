#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GraphPad Prism-style Michaelis-Menten Fitting and Calibration Curve Generator
Raw data만 입력받아 MM fitting 후 calibration curve 생성
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 백엔드 설정 (GUI 없이 PNG 저장)
import warnings
warnings.filterwarnings('ignore')


def read_raw_data(filename='prep_data/prep_raw.csv'):
    """
    prep_raw.csv에서 원본 데이터 읽기 및 정리
    
    새로운 형식:
    - 첫 번째 행: 농도 값들 (각 농도가 mean, SD, N으로 3번 반복)
    - 두 번째 행: 컬럼 헤더 (time_min, mean, SD, N, mean, SD, N, ...)
    - 세 번째 행부터: 실제 데이터
    """
    # 탭 구분자 사용
    # 첫 번째 행만 읽어서 농도 값 추출
    first_row_df = pd.read_csv(filename, header=None, nrows=1, sep='\t')
    concentration_row = first_row_df.iloc[0].values[1:]  # 첫 번째 컬럼(빈 값) 제외
    
    # 두 번째 행을 헤더로 읽기
    header_row_df = pd.read_csv(filename, header=None, skiprows=[0], nrows=1, sep='\t')
    header_names = header_row_df.iloc[0].values
    
    # 세 번째 행부터 데이터로 읽기 (헤더 없이)
    df = pd.read_csv(filename, header=None, skiprows=[0, 1], sep='\t')
    
    # 헤더 이름 설정
    df.columns = header_names
    
    # 첫 번째 컬럼이 시간
    time_col = df.columns[0]
    times = pd.to_numeric(df[time_col].values, errors='coerce')
    
    # 농도별 데이터 추출
    data = {}
    i = 1  # 첫 번째 데이터 컬럼부터 시작
    conc_idx = 0  # 농도 인덱스
    
    while i < len(df.columns):
        # 농도 값은 첫 번째 행에서 가져옴 (mean, SD, N 중 mean 위치)
        if conc_idx < len(concentration_row):
            conc_value = float(concentration_row[conc_idx * 3])  # 각 농도의 첫 번째 값 (mean 위치)
        else:
            conc_value = float(concentration_row[conc_idx * 3]) if len(concentration_row) > conc_idx * 3 else conc_idx
        
        # 컬럼명 생성
        conc_name = f"{conc_value} ug/mL"
        
        # mean 컬럼 (값)
        value_col_idx = i
        # SD 컬럼
        sd_col_idx = i + 1 if i + 1 < len(df.columns) else None
        # N 컬럼 (사용 안 함, 건너뜀)
        
        value_col = df.columns[value_col_idx]
        sd_col = df.columns[sd_col_idx] if sd_col_idx is not None else None
        
        # NaN이 아닌 값만 추출 (0 값도 포함)
        # 컬럼 인덱스로 직접 접근하여 Series 추출
        value_series = df.iloc[:, value_col_idx]
        valid_mask = (~pd.isna(value_series)) & (value_series.astype(str) != '') & (value_series.astype(str) != 'nan')
        valid_mask = valid_mask.values.flatten() if hasattr(valid_mask, 'values') else np.array(valid_mask).flatten()
        
        # 유효한 행만 필터링 (인덱스로 직접 접근)
        valid_indices = np.where(valid_mask)[0]
        valid_times = pd.to_numeric(df.iloc[valid_indices, 0], errors='coerce').values
        valid_values = pd.to_numeric(df.iloc[valid_indices, value_col_idx], errors='coerce').values
        if sd_col_idx is not None:
            valid_sd = pd.to_numeric(df.iloc[valid_indices, sd_col_idx], errors='coerce').values
        else:
            valid_sd = None
        
        # NaN만 제거 (0 값은 유지)
        valid_mask2 = ~pd.isna(valid_values)
        valid_mask2 = np.array(valid_mask2)  # numpy array로 변환
        valid_times = valid_times[valid_mask2]
        valid_values = valid_values[valid_mask2]
        if valid_sd is not None:
            valid_sd = valid_sd[valid_mask2]
        
        if len(valid_times) > 0:
            data[conc_name] = {
                'time': valid_times,
                'value': valid_values,
                'SD': valid_sd,
                'concentration': conc_value,
                'conc_name': conc_name
            }
        
        # 다음 농도로 (3개 컬럼씩: mean, SD, N)
        i += 3
        conc_idx += 1
    
    return data


def exponential_association(t, F0, Fmax, k):
    """
    Exponential Association 모델 (GraphPad Prism 표준)
    F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
    """
    return F0 + (Fmax - F0) * (1 - np.exp(-k * t))


def michaelis_menten_kinetic(t, Vmax, Km, F0):
    """
    Michaelis-Menten Kinetic 모델 (시간 도메인)
    v = Vmax * S / (Km + S)
    적분형: F(t) = F0 + Vmax * t * S / (Km + S)
    
    단순화: F(t) = F0 + (Vmax * t) / (1 + Km/t)
    또는 더 정확히는 수치 적분 필요
    """
    # 근사: F(t) = F0 + Vmax * t / (1 + Km_eff / t)
    # Km_eff는 Km을 시간 단위로 변환
    if Km > 0:
        # 더 간단한 근사 사용
        rate = Vmax / (1 + Km / t) if t > 0 else 0
        return F0 + rate * t
    else:
        return F0 + Vmax * t


def fit_time_course(times, values, model='exponential'):
    """
    시간 경과 곡선에 모델 피팅
    
    Parameters:
    - times: 시간 배열
    - values: 형광값 배열
    - model: 'exponential' 또는 'mm_kinetic'
    
    Returns:
    - params: 피팅 파라미터 딕셔너리
    - fit_values: 피팅된 값
    - r_squared: 결정계수
    """
    times = np.array(times)
    values = np.array(values)
    
    # 초기값 추정
    F0_init = values[0] if len(values) > 0 else 0
    Fmax_init = np.max(values)
    k_init = 0.1  # 초기 추정
    
    if model == 'exponential':
        # Exponential Association 모델
        try:
            # 매우 넓은 bounds 사용
            popt, pcov = curve_fit(
                exponential_association, times, values,
                p0=[F0_init, Fmax_init, k_init],
                bounds=([-1000, F0_init, 0.001], [Fmax_init, Fmax_init * 3, 10]),
                maxfev=5000
            )
            F0, Fmax, k = popt
            
            # Vmax와 Km으로 변환 (GraphPad Prism 스타일)
            # Vmax는 초기 속도, Km은 반속도 관련
            Vmax = k * (Fmax - F0)  # 초기 속도
            Km = (Fmax - F0) / 2  # 반속도 지점 근사
            
            fit_values = exponential_association(times, F0, Fmax, k)
            
        except Exception as e:
            print(f"   ⚠️ 피팅 실패: {e}, 기본값 사용")
            F0, Fmax, k = F0_init, Fmax_init, k_init
            Vmax = k * (Fmax - F0)
            Km = (Fmax - F0) / 2
            fit_values = values
    
    else:  # mm_kinetic
        try:
            popt, pcov = curve_fit(
                michaelis_menten_kinetic, times, values,
                p0=[(Fmax_init - F0_init) / times[-1], 1.0, F0_init],
                bounds=([0, 0.01, 0], [np.inf, 100, F0_init * 2]),
                maxfev=5000
            )
            Vmax, Km, F0 = popt
            Fmax = np.max(values)
            k = Vmax / (Fmax - F0) if (Fmax - F0) > 0 else 0.1
            fit_values = michaelis_menten_kinetic(times, Vmax, Km, F0)
        except:
            F0, Fmax, k = F0_init, Fmax_init, k_init
            Vmax = k * (Fmax - F0)
            Km = (Fmax - F0) / 2
            fit_values = values
    
    # R² 계산
    ss_res = np.sum((values - fit_values) ** 2)
    ss_tot = np.sum((values - np.mean(values)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    params = {
        'F0': F0,
        'Fmax': Fmax,
        'k': k,
        'Vmax': Vmax,
        'Km': Km,
        'R_squared': r_squared
    }
    
    return params, fit_values, r_squared


def michaelis_menten_calibration(x, Vmax_cal, Km_cal):
    """
    Calibration Curve: Michaelis-Menten 방정식
    y = (Vmax * x) / (Km + x)
    
    Parameters:
    - x: 농도
    - Vmax_cal: 최대 응답
    - Km_cal: 반속도 농도 (Michaelis 상수)
    """
    return (Vmax_cal * x) / (Km_cal + x)


def fit_calibration_curve(concentrations, responses):
    """
    농도 vs 응답 데이터에 MM calibration curve 피팅
    
    Parameters:
    - concentrations: 농도 배열
    - responses: 응답 배열 (Vmax 또는 형광값)
    
    Returns:
    - params: Vmax_cal, Km_cal
    - fit_values: 피팅된 값
    - equation: 방정식 문자열
    """
    concentrations = np.array(concentrations)
    responses = np.array(responses)
    
    # 초기값 추정
    Vmax_init = np.max(responses)
    Km_init = np.mean(concentrations)
    
    try:
        popt, pcov = curve_fit(
            michaelis_menten_calibration,
            concentrations, responses,
            p0=[Vmax_init, Km_init],
            bounds=([0, 0.01], [np.inf, np.inf]),
            maxfev=5000
        )
        
        Vmax_cal, Km_cal = popt
        perr = np.sqrt(np.diag(pcov))
        
        fit_values = michaelis_menten_calibration(concentrations, Vmax_cal, Km_cal)
        
        # R² 계산
        ss_res = np.sum((responses - fit_values) ** 2)
        ss_tot = np.sum((responses - np.mean(responses)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        equation = f"y = ({Vmax_cal:.2f} * x) / ({Km_cal:.4f} + x)"
        
        params = {
            'Vmax_cal': Vmax_cal,
            'Km_cal': Km_cal,
            'Vmax_cal_std': perr[0],
            'Km_cal_std': perr[1],
            'R_squared': r_squared
        }
        
        return params, fit_values, equation
        
    except Exception as e:
        print(f"   ⚠️ Calibration curve 피팅 실패: {e}")
        # 선형 근사
        coeffs = np.polyfit(concentrations, responses, 1)
        fit_values = np.polyval(coeffs, concentrations)
        equation = f"y = {coeffs[0]:.2f} * x + {coeffs[1]:.2f}"
        
        params = {
            'Vmax_cal': coeffs[0],
            'Km_cal': 0,
            'R_squared': 0
        }
        
        return params, fit_values, equation


def main():
    """메인 함수"""
    print("📊 GraphPad Prism-style MM Fitting & Calibration Curve Generator")
    print("=" * 70)
    
    # 1. Raw data 읽기
    print("\n1️⃣ Raw data 파일 읽는 중...")
    try:
        raw_data = read_raw_data('prep_data/prep_raw.csv')
        print(f"   ✅ {len(raw_data)}개 농도 조건 발견")
        for conc_name, data in raw_data.items():
            print(f"      - {conc_name}: {len(data['time'])}개 데이터 포인트")
    except Exception as e:
        print(f"   ❌ 오류: {e}")
        return
    
    # 2. 각 농도별 시간 경과 곡선 피팅
    print("\n2️⃣ 각 농도별 시간 경과 곡선 피팅 (MM/Exponential)...")
    
    mm_results = {}
    all_fit_data = []
    
    for conc_name, data in raw_data.items():
        times = data['time']
        values = data['value']
        
        # Exponential Association 모델로 피팅
        params, fit_values, r_sq = fit_time_course(times, values, model='exponential')
        
        # MM 파라미터 추출
        Vmax = params['Vmax']
        Km = params['Km']
        F0 = params['F0']
        Fmax = params['Fmax']
        
        mm_results[conc_name] = {
            'concentration': data['concentration'],
            'Vmax': Vmax,
            'Km': Km,
            'F0': F0,
            'Fmax': Fmax,
            'k': params['k'],
            'R_squared': r_sq
        }
        
        # Fit curve 데이터 저장
        for t, val, fit_val in zip(times, values, fit_values):
            all_fit_data.append({
                'Concentration': conc_name,
                'Conc_Value': data['concentration'],
                'Time_min': t,
                'Observed_Value': val,
                'Fit_Value': fit_val,
                'Residual': val - fit_val
            })
        
        print(f"   ✅ {conc_name}: Vmax={Vmax:.2f}, Km={Km:.4f}, R²={r_sq:.4f}")
    
    # 3. MM Results CSV 저장
    print("\n3️⃣ MM Results CSV 생성 중...")
    
    results_data = []
    for conc_name, params in sorted(mm_results.items(), key=lambda x: x[1]['concentration']):
        results_data.append({
            'Concentration': conc_name,
            'Conc_Value': params['concentration'],
            'Vmax': params['Vmax'],
            'Km': params['Km'],
            'F0': params['F0'],
            'Fmax': params['Fmax'],
            'k': params['k'],
            'R_squared': params['R_squared']
        })
    
    results_df = pd.DataFrame(results_data)
    results_filename = 'prep_data/MM_results_generated.csv'
    
    # GraphPad Prism 스타일로 저장
    with open(results_filename, 'w', newline='', encoding='utf-8-sig') as f:
        f.write(',')
        f.write(','.join(results_df['Concentration'].astype(str)) + '\n')
        f.write('Michaelis-Menten,\n')
        f.write('Best-fit values,\n')
        f.write(f"Vmax,{','.join(results_df['Vmax'].astype(str).str[:10])}\n")
        f.write(f"Km,{','.join(results_df['Km'].astype(str))}\n")
    
    results_df.to_csv('prep_data/MM_results_detailed.csv', index=False)
    print(f"   ✅ {results_filename} 저장 완료")
    print(f"   ✅ prep_data/MM_results_detailed.csv 저장 완료 (상세 데이터)")
    
    # 4. Calibration Curve 생성
    print("\n4️⃣ Calibration Curve 생성 중...")
    
    # 농도 vs Vmax로 calibration curve 피팅
    concentrations = [mm_results[cn]['concentration'] for cn in sorted(mm_results.keys(), 
                                                                      key=lambda x: mm_results[x]['concentration'])]
    vmax_values = [mm_results[cn]['Vmax'] for cn in sorted(mm_results.keys(), 
                                                          key=lambda x: mm_results[x]['concentration'])]
    
    # MM calibration curve 피팅
    cal_params, cal_fit_values, cal_equation = fit_calibration_curve(concentrations, vmax_values)
    
    print(f"   ✅ Calibration Equation: {cal_equation}")
    print(f"      Vmax_cal = {cal_params['Vmax_cal']:.2f} ± {cal_params.get('Vmax_cal_std', 0):.2f}")
    print(f"      Km_cal = {cal_params['Km_cal']:.4f} ± {cal_params.get('Km_cal_std', 0):.4f}")
    print(f"      R² = {cal_params['R_squared']:.4f}")
    
    # 5. Calibration Curve XY 데이터 생성
    print("\n5️⃣ Calibration Curve XY 데이터 생성 중...")
    
    # 고밀도 농도 범위
    conc_min = min(concentrations)
    conc_max = max(concentrations)
    conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
    
    # Calibration curve 계산
    cal_y_values = michaelis_menten_calibration(conc_range, 
                                                cal_params['Vmax_cal'], 
                                                cal_params['Km_cal'])
    
    # Calibration curve 데이터 저장
    cal_curve_data = []
    for x, y in zip(conc_range, cal_y_values):
        cal_curve_data.append({
            'Concentration': x,
            'Vmax_Fitted': y,
            'Equation': cal_equation
        })
    
    cal_curve_df = pd.DataFrame(cal_curve_data)
    cal_curve_filename = 'prep_data/MM_calibration_curve.csv'
    cal_curve_df.to_csv(cal_curve_filename, index=False)
    print(f"   ✅ {cal_curve_filename} 저장 완료 ({len(cal_curve_df)} 행)")
    
    # 6. Fit curves 데이터 저장
    fit_curves_df = pd.DataFrame(all_fit_data)
    fit_curves_filename = 'prep_data/MM_fit_curves.csv'
    fit_curves_df.to_csv(fit_curves_filename, index=False)
    print(f"   ✅ {fit_curves_filename} 저장 완료 ({len(fit_curves_df)} 행)")
    
    # 7. 방정식 요약 저장
    print("\n6️⃣ 방정식 요약 저장 중...")
    
    equations_data = [{
        'Type': 'Calibration Curve',
        'Equation': cal_equation,
        'Vmax': cal_params['Vmax_cal'],
        'Km': cal_params['Km_cal'],
        'R_squared': cal_params['R_squared']
    }]
    
    # 각 농도별 시간 곡선 방정식
    for conc_name, params in sorted(mm_results.items(), key=lambda x: x[1]['concentration']):
        eq = f"F(t) = {params['F0']:.2f} + ({params['Fmax'] - params['F0']:.2f}) * [1 - exp(-{params['k']:.4f}*t)]"
        equations_data.append({
            'Type': f'Time Course ({conc_name})',
            'Equation': eq,
            'Vmax': params['Vmax'],
            'Km': params['Km'],
            'R_squared': params['R_squared']
        })
    
    equations_df = pd.DataFrame(equations_data)
    equations_filename = 'prep_data/MM_equations.csv'
    equations_df.to_csv(equations_filename, index=False)
    print(f"   ✅ {equations_filename} 저장 완료")
    
    # MM_calibration_equations.csv 형식으로 저장 (사용자 요청 형식)
    calibration_equations_data = []
    for conc_name, params in sorted(mm_results.items(), key=lambda x: x[1]['concentration']):
        # F0가 거의 0이면 단순화된 형식 사용
        if abs(params['F0']) < 0.01:
            eq_str = f"F(t) = 0.00 + ({params['Fmax']:.2f}) * [1 - exp(-{params['k']:.4f}*t)]"
        else:
            eq_str = f"F(t) = {params['F0']:.2f} + ({params['Fmax']:.2f}) * [1 - exp(-{params['k']:.4f}*t)]"
        
        calibration_equations_data.append({
            'Concentration': conc_name,
            'F0': params['F0'],
            'Fmax': params['Fmax'],
            'Vmax': params['Vmax'],
            'Km': params['Km'],
            'k_rate': params['k'],
            'Equation': eq_str
        })
    
    calibration_equations_df = pd.DataFrame(calibration_equations_data)
    calibration_equations_filename = 'prep_data/MM_calibration_equations.csv'
    calibration_equations_df.to_csv(calibration_equations_filename, index=False)
    print(f"   ✅ {calibration_equations_filename} 저장 완료 (농도별 상세 방정식)")
    
    # 최종 요약
    print("\n" + "=" * 70)
    print("📋 생성된 파일:")
    print(f"   1. {results_filename} - GraphPad Prism 스타일 MM 결과")
    print(f"   2. prep_data/MM_results_detailed.csv - 상세 MM 파라미터")
    print(f"   3. {cal_curve_filename} - Calibration curve XY 데이터 (그래프용)")
    print(f"   4. prep_data/MM_calibration_curve.png - Calibration curve 그래프 (PNG)")
    print(f"   5. {fit_curves_filename} - 각 농도별 시간 곡선 fit 데이터")
    print(f"   6. {equations_filename} - 모든 방정식 요약")
    print("\n📊 Calibration Curve:")
    print(f"   {cal_equation}")
    print(f"   농도 범위: {conc_min:.4f} - {conc_max:.4f} (확장: {conc_min*0.5:.4f} - {conc_max*1.5:.4f})")
    # 7. Calibration Curve 그래프 생성
    print("\n7️⃣ Calibration Curve 그래프 생성 중...")
    plot_calibration_curve(
        cal_curve_df, results_df, cal_params, cal_equation
    )
    print("   ✅ prep_data/MM_calibration_curve.png 저장 완료")
    
    print("\n✨ 완료!")


def plot_calibration_curve(cal_curve_df, results_df, cal_params, cal_equation):
    """
    Calibration curve 그래프를 그리고 PNG로 저장
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Calibration curve 그리기
    ax.plot(
        cal_curve_df['Concentration'],
        cal_curve_df['Vmax_Fitted'],
        'b-', linewidth=2.5,
        label=f'MM Fit: {cal_equation}',
        zorder=1
    )
    
    # 실험 데이터 포인트 그리기
    concentrations = results_df['Conc_Value'].values
    vmax_values = results_df['Vmax'].values
    
    ax.scatter(
        concentrations,
        vmax_values,
        color='red',
        s=150,
        marker='o',
        edgecolors='black',
        linewidths=2,
        label='Experimental Data',
        zorder=2
    )
    
    # 그래프 스타일
    ax.set_xlabel('Concentration (μg/mL)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Vmax (Fluorescence Units)', fontsize=14, fontweight='bold')
    ax.set_title('Michaelis-Menten Calibration Curve', fontsize=16, fontweight='bold', pad=20)
    
    # 그리드 추가
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=12, loc='lower right', framealpha=0.9)
    
    # 통계 정보 텍스트 박스
    stats_text = f"Vmax_cal = {cal_params['Vmax_cal']:.2f}\n"
    stats_text += f"Km_cal = {cal_params['Km_cal']:.4f}\n"
    stats_text += f"R² = {cal_params['R_squared']:.4f}"
    
    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # 레이아웃 조정
    plt.tight_layout()
    
    # PNG 저장
    plt.savefig('prep_data/MM_calibration_curve.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()
