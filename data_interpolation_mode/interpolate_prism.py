#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interpolation for RFU vs Time curves
Prism 보간 알고리즘 구현
"""

import pandas as pd
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import warnings
warnings.filterwarnings('ignore')


def exponential_association(t, F0, Fmax, k):
    """
    Exponential Association 모델
    F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
    """
    return F0 + (Fmax - F0) * (1 - np.exp(-k * t))


def inverse_exponential_association(y, F0, Fmax, k, x_low, x_high):
    """
    Y → X 역함수 (이진 이등분법 사용)
    F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
    → t = -ln(1 - (y - F0)/(Fmax - F0)) / k
    
    단, 단조성이 보장되지 않는 경우를 대비해 이진 이등분법 사용
    """
    def f(x):
        return exponential_association(x, F0, Fmax, k) - y
    
    try:
        # 단조 함수인 경우 직접 계산
        if Fmax > F0 and k > 0:
            if y <= F0:
                return x_low
            if y >= Fmax:
                return x_high
            
            # 단조 구간 내에서 직접 계산 시도
            t_direct = -np.log(1 - (y - F0) / (Fmax - F0)) / k
            if x_low <= t_direct <= x_high:
                return t_direct
        
        # 이진 이등분법 (brentq 사용)
        return brentq(f, x_low, x_high, xtol=1e-10, maxiter=100)
    except:
        # 실패 시 선형 보간
        y_low = exponential_association(x_low, F0, Fmax, k)
        y_high = exponential_association(x_high, F0, Fmax, k)
        if abs(y_high - y_low) < 1e-10:
            return x_low
        return x_low + (y - y_low) / (y_high - y_low) * (x_high - x_low)


def create_prism_interpolation_range(x_data):
    """
    Prism 스타일 보간/외삽 범위 생성
    
    규칙:
    1. 데이터 범위: [Xmin, Xmax]
    2. 확장 범위: 각 방향으로 (Xmax - Xmin) / 2
    3. 모든 데이터가 양수면 음수 제외
    4. 모든 데이터가 음수면 양수 제외
    """
    x_min = np.min(x_data)
    x_max = np.max(x_data)
    x_range = x_max - x_min
    
    # 확장 거리
    extension = x_range / 2
    
    # 확장된 범위
    interp_min = x_min - extension
    interp_max = x_max + extension
    
    # 특별한 경우 처리
    if np.all(x_data >= 0):
        interp_min = max(0, interp_min)
    
    if np.all(x_data <= 0):
        interp_max = min(0, interp_max)
    
    return interp_min, interp_max


def prism_interpolate_x_to_y(x_interp, F0, Fmax, k):
    """
    X → Y 보간 (방정식으로 직접 계산)
    """
    return exponential_association(x_interp, F0, Fmax, k)


def prism_interpolate_y_to_x(y_interp, F0, Fmax, k, x_range_min, x_range_max, 
                               x_data_min, x_data_max, n_segments=1000):
    """
    Y → X 보간 (Prism 알고리즘)
    
    1. 확장된 범위를 1000개 선분으로 나눔
    2. 가장 낮은 X부터 스캔하여 Y가 포함된 첫 번째 선분 찾기
    3. 해당 선분 내에서 이진 이등분법으로 정확한 X 계산
    """
    # 확장된 범위를 선분으로 나눔
    segment_edges = np.linspace(x_range_min, x_range_max, n_segments + 1)
    
    # 각 선분의 Y 값 계산
    segment_y_low = exponential_association(segment_edges[:-1], F0, Fmax, k)
    segment_y_high = exponential_association(segment_edges[1:], F0, Fmax, k)
    
    # Y가 포함된 선분 찾기 (가장 낮은 X부터)
    x_result = []
    
    for y in y_interp:
        found = False
        
        # 데이터 범위 내에서 먼저 찾기
        for i in range(len(segment_edges) - 1):
            seg_low = segment_edges[i]
            seg_high = segment_edges[i + 1]
            
            # 선분이 데이터 범위 내에 있는지 확인
            if seg_high < x_data_min or seg_low > x_data_max:
                continue
            
            y_low = segment_y_low[i]
            y_high = segment_y_high[i]
            
            # Y가 이 선분에 포함되는지 확인
            if (min(y_low, y_high) <= y <= max(y_low, y_high)):
                # 이진 이등분법으로 정확한 X 계산
                try:
                    x_found = inverse_exponential_association(
                        y, F0, Fmax, k, seg_low, seg_high
                    )
                    x_result.append(x_found)
                    found = True
                    break
                except:
                    continue
        
        # 데이터 범위 내에서 못 찾으면 외삽 시도
        if not found:
            # X < Xmin 범위
            for i in range(len(segment_edges) - 1):
                if segment_edges[i + 1] >= x_data_min:
                    break
                
                seg_low = segment_edges[i]
                seg_high = segment_edges[i + 1]
                y_low = segment_y_low[i]
                y_high = segment_y_high[i]
                
                if (min(y_low, y_high) <= y <= max(y_low, y_high)):
                    try:
                        x_found = inverse_exponential_association(
                            y, F0, Fmax, k, seg_low, seg_high
                        )
                        x_result.append(x_found)
                        found = True
                        break
                    except:
                        continue
            
            # X > Xmax 범위
            if not found:
                for i in range(len(segment_edges) - 2, -1, -1):
                    if segment_edges[i] <= x_data_max:
                        break
                    
                    seg_low = segment_edges[i]
                    seg_high = segment_edges[i + 1]
                    y_low = segment_y_low[i]
                    y_high = segment_y_high[i]
                    
                    if (min(y_low, y_high) <= y <= max(y_low, y_high)):
                        try:
                            x_found = inverse_exponential_association(
                                y, F0, Fmax, k, seg_low, seg_high
                            )
                            x_result.append(x_found)
                            found = True
                            break
                        except:
                            continue
        
        # 여전히 못 찾으면 NaN
        if not found:
            x_result.append(np.nan)
    
    return np.array(x_result)


def main():
    """메인 함수"""
    print("📊 GraphPad Interpolation Generator")
    print("=" * 70)
    
    # 1. Fit 파라미터 읽기
    print("\n1️⃣ Fit 파라미터 읽는 중...")
    try:
        mm_results = pd.read_csv('MM_results_detailed.csv')
        print(f"   ✅ {len(mm_results)}개 농도 조건 발견")
    except Exception as e:
        print(f"   ❌ 오류: {e}")
        print("   💡 먼저 prep.py를 실행하여 MM_results_detailed.csv를 생성해주세요.")
        return
    
    # 2. Raw data 읽기
    print("\n2️⃣ Raw data 읽는 중...")
    try:
        raw_df = pd.read_csv('prep_raw.csv', header=0)
        time_col = raw_df.columns[0]
        times = raw_df[time_col].values
        print(f"   ✅ 시간 범위: {times.min():.1f} - {times.max():.1f} 분")
    except Exception as e:
        print(f"   ❌ 오류: {e}")
        return
    
    # 3. 각 농도별로 보간 데이터 생성
    print("\n3️⃣ Prism 스타일 보간 데이터 생성 중...")
    
    all_interp_data = []
    x_interp_all = []
    y_interp_all = []
    
    # 보간 범위 결정 (시간 축)
    x_data_min = times.min()
    x_data_max = times.max()
    x_range_min, x_range_max = create_prism_interpolation_range(times)
    
    print(f"   보간 범위: {x_range_min:.3f} - {x_range_max:.3f} 분")
    print(f"   데이터 범위: {x_data_min:.1f} - {x_data_max:.1f} 분")
    
    # 고밀도 X 포인트 생성 (1000개 선분)
    n_points = 1000
    x_interp = np.linspace(x_range_min, x_range_max, n_points + 1)
    x_interp_all = x_interp
    
    for idx, row in mm_results.iterrows():
        conc_name = row['Concentration']
        F0 = row['F0']
        Fmax = row['Fmax']
        k = row['k']
        
        print(f"\n   📊 {conc_name}:")
        print(f"      F0={F0:.2f}, Fmax={Fmax:.2f}, k={k:.4f}")
        
        # X → Y 보간 (고밀도 곡선)
        y_interp = prism_interpolate_x_to_y(x_interp, F0, Fmax, k)
        y_interp_all.append(y_interp)
        
        # 보간 데이터 저장
        for x, y in zip(x_interp, y_interp):
            all_interp_data.append({
                'Concentration': conc_name,
                'Concentration [ug/mL]': row['Concentration [ug/mL]'],
                'Time_min': x,
                'RFU_Interpolated': y,
                'Is_Extrapolated': (x < x_data_min) or (x > x_data_max)
            })
        
        print(f"      ✅ {len(x_interp)}개 포인트 생성 완료")
    
    # 4. CSV 저장
    print("\n4️⃣ 보간 데이터 저장 중...")
    
    interp_df = pd.DataFrame(all_interp_data)
    interp_filename = 'MM_interpolated_curves.csv'
    interp_df.to_csv(interp_filename, index=False)
    print(f"   ✅ {interp_filename} 저장 완료 ({len(interp_df)} 행)")
    
    # 5. 요약 테이블 생성
    print("\n5️⃣ 요약 테이블 생성 중...")
    
    summary_data = []
    for idx, row in mm_results.iterrows():
        conc_name = row['Concentration']
        subset = interp_df[interp_df['Concentration'] == conc_name]
        
        interp_range = subset[
            (subset['Time_min'] >= x_data_min) & 
            (subset['Time_min'] <= x_data_max)
        ]
        extrap_range = subset[
            (subset['Time_min'] < x_data_min) | 
            (subset['Time_min'] > x_data_max)
        ]
        
        summary_data.append({
            'Concentration': conc_name,
            'Interpolation_Range_min': x_range_min,
            'Interpolation_Range_max': x_range_max,
            'Data_Range_min': x_data_min,
            'Data_Range_max': x_data_max,
            'Num_Interpolation_Points': len(interp_range),
            'Num_Extrapolation_Points': len(extrap_range),
            'Total_Points': len(subset)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_filename = 'MM_interpolation_summary.csv'
    summary_df.to_csv(summary_filename, index=False)
    print(f"   ✅ {summary_filename} 저장 완료")
    
    # 6. Y → X 보간 예제 (선택적)
    print("\n6️⃣ Y → X 보간 예제 생성 중...")
    
    # 각 농도별로 몇 가지 Y 값에 대해 X 계산
    y_to_x_examples = []
    
    for idx, row in mm_results.iterrows():
        conc_name = row['Concentration']
        F0 = row['F0']
        Fmax = row['Fmax']
        k = row['k']
        
        # Y 값 예제 (F0에서 Fmax까지 몇 개)
        y_examples = np.linspace(F0 + (Fmax - F0) * 0.1, 
                                Fmax - (Fmax - F0) * 0.1, 5)
        
        x_calculated = prism_interpolate_y_to_x(
            y_examples, F0, Fmax, k,
            x_range_min, x_range_max,
            x_data_min, x_data_max
        )
        
        for y, x in zip(y_examples, x_calculated):
            if not np.isnan(x):
                y_to_x_examples.append({
                    'Concentration': conc_name,
                    'Target_RFU': y,
                    'Calculated_Time_min': x,
                    'Is_In_Data_Range': (x_data_min <= x <= x_data_max)
                })
    
    if y_to_x_examples:
        y_to_x_df = pd.DataFrame(y_to_x_examples)
        y_to_x_filename = 'MM_Y_to_X_interpolation.csv'
        y_to_x_df.to_csv(y_to_x_filename, index=False)
        print(f"   ✅ {y_to_x_filename} 저장 완료 ({len(y_to_x_df)} 행)")
    
    # 7. 보간 곡선 그래프 생성
    print("\n7️⃣ 보간 곡선 그래프 생성 중...")
    plot_interpolated_curves(interp_df, raw_df, mm_results)
    print("   ✅ MM_interpolated_curves.png 저장 완료")
    
    # 최종 요약
    print("\n" + "=" * 70)
    print("📋 생성된 파일:")
    print(f"   1. {interp_filename} - Prism 스타일 보간 곡선 데이터")
    print(f"   2. {summary_filename} - 보간 범위 요약")
    if y_to_x_examples:
        print(f"   3. {y_to_x_filename} - Y → X 보간 예제")
    print(f"   4. MM_interpolated_curves.png - 보간 곡선 그래프")
    print("\n📊 보간 정보:")
    print(f"   보간 범위: {x_range_min:.3f} - {x_range_max:.3f} 분")
    print(f"   선분 개수: {n_points}")
    print(f"   데이터 범위: {x_data_min:.1f} - {x_data_max:.1f} 분")
    print("\n✨ 완료!")


def plot_interpolated_curves(interp_df, raw_df, mm_results):
    """
    보간 곡선과 원본 데이터를 함께 그래프로 그리기
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 농도별 색상 매핑 (농도 순서대로)
    conc_color_map = {
        '0.3125ug/ml': 'blue',
        '0.625 ug/mL': 'red',
        '1.25 ug/mL': 'orange',
        '2.5 ug/mL': 'green',
        '5 ug/mL': 'purple'
    }
    
    # 농도별 컬럼 인덱스 매핑 (prep_raw.csv 구조: 시간, 값, SD, 복제수)
    conc_col_map = {
        '0.3125ug/ml': 1,
        '0.625 ug/mL': 4,
        '1.25 ug/mL': 7,
        '2.5 ug/mL': 10,
        '5 ug/mL': 13
    }
    
    # 농도 순서 정렬
    conc_order = ['0.3125ug/ml', '0.625 ug/mL', '1.25 ug/mL', '2.5 ug/mL', '5 ug/mL']
    
    # 원본 데이터 읽기
    time_col = raw_df.columns[0]
    times = raw_df.iloc[1:, 0].values  # 첫 번째 행(헤더) 제외
    
    # 각 농도별로 그래프 그리기
    for conc_name in conc_order:
        if conc_name not in conc_color_map:
            continue
            
        color = conc_color_map[conc_name]
        
        # 보간 곡선
        subset = interp_df[interp_df['Concentration'] == conc_name]
        if len(subset) == 0:
            continue
            
        interp_in_range = subset[~subset['Is_Extrapolated']]
        interp_extrap = subset[subset['Is_Extrapolated']]
        
        # 데이터 범위 내 보간 곡선 (실선)
        if len(interp_in_range) > 0:
            ax.plot(interp_in_range['Time_min'], interp_in_range['RFU_Interpolated'],
                   color=color, linewidth=2.5, label=f'{conc_name} (Interpolated)',
                   zorder=2)
        
        # 외삽 영역 (점선)
        if len(interp_extrap) > 0:
            ax.plot(interp_extrap['Time_min'], interp_extrap['RFU_Interpolated'],
                   color=color, linewidth=2, linestyle='--', alpha=0.5,
                   zorder=1)
        
        # 원본 데이터 포인트 찾기 (컬럼 인덱스 직접 사용)
        if conc_name in conc_col_map:
            col_idx = conc_col_map[conc_name]
            
            if col_idx < len(raw_df.columns):
                values = raw_df.iloc[1:, col_idx].values  # 첫 번째 행(헤더) 제외
                sd_col_idx = col_idx + 1 if col_idx + 1 < len(raw_df.columns) else None
                sd_values = raw_df.iloc[1:, sd_col_idx].values if sd_col_idx else None
                
                # 데이터 포인트 플롯
                valid_mask = ~pd.isna(values) & (values > 0)
                if np.sum(valid_mask) > 0:
                    ax.scatter(times[valid_mask], values[valid_mask],
                              color=color, s=100, marker='o',
                              edgecolors='white', linewidths=1.5,
                              zorder=3, label=f'{conc_name} (Data)')
                    
                    # Error bars
                    if sd_values is not None and len(sd_values) == len(values):
                        ax.errorbar(times[valid_mask], values[valid_mask],
                                   yerr=sd_values[valid_mask],
                                   color=color, fmt='none', alpha=0.7,
                                   capsize=3, capthick=1, zorder=2)
    
    # 그래프 스타일
    ax.set_xlabel('Time (min)', fontsize=14, fontweight='bold')
    ax.set_ylabel('RFU', fontsize=14, fontweight='bold')
    ax.set_title('Time-Fluorescence Curve', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # 그리드
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # 범례
    ax.legend(fontsize=10, loc='best', framealpha=0.9, ncol=2)
    
    # 축 범위 설정
    ax.set_xlim([-2, 32])
    ax.set_ylim(bottom=0)
    
    # 레이아웃 조정
    plt.tight_layout()
    
    # PNG 저장
    plt.savefig('MM_interpolated_curves.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()

