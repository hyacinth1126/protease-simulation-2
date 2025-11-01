import os
import numpy as np
import pandas as pd
import streamlit as st
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from prep_raw_data_mode.prep import (
    read_raw_data,
    fit_time_course,
    fit_calibration_curve,
    michaelis_menten_calibration
)
from data_interpolation_mode.interpolate_prism import (
    exponential_association,
    create_prism_interpolation_range
)


def prep_raw_data_mode(st):
    """Prep Raw Data 모드 - Michaelis-Menten Fitting"""
    
    # 폴더 구조 생성
    os.makedirs("prep_data/raw", exist_ok=True)
    os.makedirs("prep_raw_data_mode/results", exist_ok=True)
    
    st.header("📊 Prep Raw Data 모드")
    st.markdown("Michaelis-Menten Fitting 및 Calculated Curve 생성")
    st.markdown("---")
    
    # 사이드바 설정
    st.sidebar.title("⚙️ Prep Raw Data 설정")
    
    # 데이터 업로드
    st.sidebar.subheader("📁 데이터 업로드")
    uploaded_file = st.sidebar.file_uploader(
        "Prep Raw CSV 파일 업로드",
        type=['csv'],
        help="prep_raw.csv 형식: 시간, 농도별 값, SD, 복제수 (3개 컬럼씩)"
    )
    
    # 샘플 데이터 다운로드
    try:
        with open("prep_data/raw/prep_raw.csv", "rb") as f:
            sample_bytes = f.read()
        st.sidebar.download_button(
            label="샘플 prep_raw.csv 다운로드",
            data=sample_bytes,
            file_name="prep_raw_sample.csv",
            mime="text/csv"
        )
    except Exception:
        pass
    
    # 데이터 로드
    if uploaded_file is not None:
        # 업로드된 파일을 임시로 저장하고 읽기
        import tempfile
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv', mode='wb') as tmp_file:
            tmp_file.write(uploaded_file.getbuffer())
            tmp_path = tmp_file.name
        
        try:
            raw_data = read_raw_data(tmp_path)
            os.unlink(tmp_path)
        except Exception as e:
            st.error(f"파일 읽기 오류: {e}")
            os.unlink(tmp_path)
            return
    else:
        # 기본 샘플 데이터 사용
        try:
            raw_data = read_raw_data('prep_data/raw/prep_raw.csv')
            st.sidebar.info("prep_data/raw/prep_raw.csv 사용 중")
        except FileNotFoundError:
            st.error("데이터 파일을 찾을 수 없습니다. CSV 파일을 업로드해주세요.")
            st.stop()
    
    # 데이터 미리보기
    st.subheader("📋 데이터 미리보기")
    
    # 반응 시간 계산 (최대값)
    all_times = [time_val for data in raw_data.values() for time_val in data['time']]
    reaction_time = f"{max(all_times):.0f} min"
    
    # N 값 읽기 (prep_raw.csv에서 직접 읽기)
    try:
        if uploaded_file is not None:
            # 업로드된 파일에서 읽기
            uploaded_file.seek(0)
            first_line = uploaded_file.readline().decode('utf-8')
            second_line = uploaded_file.readline().decode('utf-8')
            third_line = uploaded_file.readline().decode('utf-8')
            n_value = int(third_line.split('\t')[3]) if len(third_line.split('\t')) > 3 else 50
            uploaded_file.seek(0)
        else:
            # 기본 파일에서 읽기
            with open('prep_data/raw/prep_raw.csv', 'r', encoding='utf-8') as f:
                f.readline()  # 첫 번째 줄 건너뛰기
                f.readline()  # 두 번째 줄 건너뛰기
                third_line = f.readline()
                n_value = int(third_line.split('\t')[3]) if len(third_line.split('\t')) > 3 else 50
    except:
        n_value = 50  # 기본값
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("농도 조건 수", len(raw_data))
    with col2:
        st.metric("반응 시간", reaction_time)
    with col3:
        st.metric("N(시험 수)", n_value)
    
    # 농도별 정보 표시
    with st.expander("농도별 데이터 정보", expanded=False):
        # 농도 순서대로 정렬
        sorted_conc = sorted(raw_data.items(), key=lambda x: x[1]['concentration'])
        
        # 첫 번째 농도의 시간 데이터를 기준으로 사용
        first_data = sorted_conc[0][1]
        times = first_data['time']
        
        # 가로로 넓은 테이블 생성
        detail_data = {'time_min': times}
        
        for conc_name, data in sorted_conc:
            conc_label = f"{data['concentration']}"
            detail_data[f'{conc_label}_mean'] = data['value']
            
            # SD가 있으면 추가
            if data.get('SD') is not None:
                detail_data[f'{conc_label}_SD'] = data['SD']
        
        detail_df = pd.DataFrame(detail_data)
        st.dataframe(detail_df, use_container_width=True, hide_index=True, height=400)
    
    # 분석 실행 버튼
    if st.button("🚀 MM Fitting 및 Calibration Curve 생성", type="primary"):
        with st.spinner("분석 진행 중..."):
            # 진행 상황 표시
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # 1. 각 농도별 시간 경과 곡선 피팅
            status_text.text("1️⃣ 각 농도별 시간 경과 곡선 피팅 중...")
            progress_bar.progress(0.2)
            
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
                        'Concentration [ug/mL]': data['concentration'],
                        'Time_min': t,
                        'Observed_Value': val,
                        'Fit_Value': fit_val,
                        'Residual': val - fit_val
                    })
            
            progress_bar.progress(0.4)
            
            # 2. Calibration Curve 생성
            status_text.text("2️⃣ Calibration Curve 생성 중...")
            
            concentrations = [mm_results[cn]['concentration'] for cn in sorted(mm_results.keys(), 
                                                                              key=lambda x: mm_results[x]['concentration'])]
            vmax_values = [mm_results[cn]['Vmax'] for cn in sorted(mm_results.keys(), 
                                                                    key=lambda x: mm_results[x]['concentration'])]
            
            # MM calibration curve 피팅
            cal_params, cal_fit_values, cal_equation = fit_calibration_curve(concentrations, vmax_values)
            
            progress_bar.progress(0.6)
            
            # 3. 계산된 곡선 생성 (Fitting 방정식 사용)
            status_text.text("3️⃣ 계산된 곡선 생성 중...")
            
            # 시간 범위 계산
            all_times = [time_val for data in raw_data.values() for time_val in data['time']]
            x_data_min = min(all_times)
            x_data_max = max(all_times)
            x_range_min, x_range_max = create_prism_interpolation_range(np.array(all_times))
            
            # 고밀도 계산 포인트 생성 (1000개)
            n_points = 1000
            x_calc = np.linspace(x_range_min, x_range_max, n_points + 1)
            
            # 각 농도별 계산된 데이터 생성
            all_calc_data = []
            for conc_name, params in mm_results.items():
                F0 = params['F0']
                Fmax = params['Fmax']
                k = params['k']
                
                # Fitting 방정식으로 계산: F(t) = F0 + (Fmax - F0) * [1 - exp(-k*t)]
                y_calc = exponential_association(x_calc, F0, Fmax, k)
                
                for x, y in zip(x_calc, y_calc):
                    all_calc_data.append({
                        'Concentration': conc_name,
                        'Concentration [ug/mL]': params['concentration'],
                        'Time_min': x,
                        'RFU_Calculated': y,
                        'Is_Extrapolated': (x < x_data_min) or (x > x_data_max)
                    })
            
            calc_df = pd.DataFrame(all_calc_data)
            
            progress_bar.progress(0.7)
            
            # 4. 결과 데이터 준비
            status_text.text("4️⃣ 결과 준비 중...")
            
            # 결과 데이터프레임 생성 (Calibration Curve를 첫 행에 추가)
            results_data = []
            
            # Calibration Curve 추가
            results_data.append({
                'Type': 'Calibration Curve',
                'Equation': cal_equation,
                'Concentration': 'Calibration Curve',
                'Concentration [ug/mL]': None,
                'Vmax': cal_params['Vmax_cal'],
                'Km': cal_params['Km_cal'],
                'F0': None,
                'Fmax': None,
                'k': None,
                'R_squared': cal_params['R_squared']
            })
            
            # 각 농도별 Time Course 추가
            for conc_name, params in sorted(mm_results.items(), key=lambda x: x[1]['concentration']):
                eq = f"F(t) = {params['F0']:.2f} + ({params['Fmax'] - params['F0']:.2f}) * [1 - exp(-{params['k']:.4f}*t)]"
                results_data.append({
                    'Type': f'{conc_name}',
                    'Equation': eq,
                    'Concentration': conc_name,
                    'Concentration [ug/mL]': params['concentration'],
                    'Vmax': params['Vmax'],
                    'Km': params['Km'],
                    'F0': params['F0'],
                    'Fmax': params['Fmax'],
                    'k': params['k'],
                    'R_squared': params['R_squared']
                })
            
            results_df = pd.DataFrame(results_data)
            fit_curves_df = pd.DataFrame(all_fit_data)
            
            # Calibration curve 데이터
            conc_min = min(concentrations)
            conc_max = max(concentrations)
            conc_range = np.linspace(conc_min * 0.5, conc_max * 1.5, 200)
            cal_y_values = michaelis_menten_calibration(conc_range, 
                                                        cal_params['Vmax_cal'], 
                                                        cal_params['Km_cal'])
            
            cal_curve_df = pd.DataFrame({
                'Concentration_ug/mL': conc_range,
                'Vmax_Fitted': cal_y_values,
                'Equation': cal_equation
            })
            
            progress_bar.progress(0.95)
            status_text.text("💾 결과 파일 저장 중...")
            
            # 결과 파일을 fitting_results 폴더에 자동 저장
            try:
                # Calculated curves 저장
                calc_df.to_csv('prep_raw_data_mode/results/MM_calculated_curves.csv', index=False)
                
                # MM results 저장
                results_df.to_csv('prep_raw_data_mode/results/MM_results_detailed.csv', index=False)
                
                # Fit curves 저장
                fit_curves_df.to_csv('prep_raw_data_mode/results/MM_fit_curves.csv', index=False)
                
                # Calibration curve 저장
                cal_curve_df.to_csv('prep_raw_data_mode/results/MM_calibration_curve.csv', index=False)
                
                st.sidebar.success("✅ 결과 파일이 prep_raw_data_mode/results/ 에 저장되었습니다!")
            except Exception as e:
                st.sidebar.warning(f"⚠️ 파일 저장 중 오류: {e}")
            
            progress_bar.progress(1.0)
            status_text.text("✅ 분석 완료!")
            
            # Session state에 저장
            st.session_state['prep_results'] = {
                'mm_results': mm_results,
                'results_df': results_df,
                'fit_curves_df': fit_curves_df,
                'calc_df': calc_df,
                'x_data_min': x_data_min,
                'x_data_max': x_data_max,
                'cal_params': cal_params,
                'cal_equation': cal_equation,
                'cal_curve_df': cal_curve_df,
                'raw_data': raw_data
            }
    
    # 결과 표시
    if 'prep_results' in st.session_state:
        results = st.session_state['prep_results']
        
        # 탭 구성
        tab1, tab2, tab3 = st.tabs([
            "📈 Michaelis-Menten Calibration Curve",
            "📊 Michaelis-Menten Calibration Results",
            "💾 Download"
        ])
        
        with tab1:
            st.caption("F(t) = F₀ + (Fₘₐₓ - F₀) × [1 - exp(-k·t)]")
            
            # 각 농도별 그래프
            fig = make_subplots(
                rows=1, cols=1,
                subplot_titles=('Time-Fluorescence Curve',)
            )
            
            colors = ['blue', 'red', 'orange', 'green', 'purple']
            conc_order = sorted(results['results_df']['Concentration'].values, 
                              key=lambda x: results['results_df'][results['results_df']['Concentration']==x]['Concentration [ug/mL]'].values[0])
            
            x_data_min = results['x_data_min']
            x_data_max = results['x_data_max']
            
            for idx, conc_name in enumerate(conc_order):
                color = colors[idx % len(colors)]
                
                # 계산된 곡선 (Fitting 방정식 사용)
                calc_subset = results['calc_df'][results['calc_df']['Concentration'] == conc_name]
                
                # 데이터 범위 내 곡선 (실선)
                calc_in_range = calc_subset[~calc_subset['Is_Extrapolated']]
                if len(calc_in_range) > 0:
                    fig.add_trace(go.Scatter(
                        x=calc_in_range['Time_min'],
                        y=calc_in_range['RFU_Calculated'],
                        mode='lines',
                        name=f'{conc_name} (Fitted)',
                        line=dict(color=color, width=2.5),
                        legendgroup=conc_name,
                        showlegend=True
                    ))
                
                # 외삽 영역 (점선)
                calc_extrap = calc_subset[calc_subset['Is_Extrapolated']]
                if len(calc_extrap) > 0:
                    fig.add_trace(go.Scatter(
                        x=calc_extrap['Time_min'],
                        y=calc_extrap['RFU_Calculated'],
                        mode='lines',
                        name=f'{conc_name} (Extrapolated)',
                        line=dict(color=color, width=2, dash='dash'),
                        opacity=0.5,
                        legendgroup=conc_name,
                        showlegend=False
                    ))
                
                # 원본 데이터 포인트
                raw_subset = results['fit_curves_df'][results['fit_curves_df']['Concentration'] == conc_name]
                if len(raw_subset) > 0:
                    fig.add_trace(go.Scatter(
                        x=raw_subset['Time_min'],
                        y=raw_subset['Observed_Value'],
                        mode='markers',
                        name=f'{conc_name} (Data)',
                        marker=dict(color=color, size=10, 
                                   line=dict(color='white', width=1.5)),
                        legendgroup=conc_name,
                        showlegend=True
                    ))
            
            fig.update_layout(
                xaxis_title='Time (min)',
                yaxis_title='RFU',
                height=700,
                template='plotly_white',
                hovermode='x unified',
                legend=dict(
                    orientation="v",
                    yanchor="top",
                    y=0.99,
                    xanchor="right",
                    x=0.99,
                    bgcolor="rgba(255,255,255,0.8)"
                )
            )
            
            fig.update_xaxes(range=[-2, x_data_max + 2])
            fig.update_yaxes(rangemode='tozero')
            
            st.plotly_chart(fig, use_container_width=True)
            
            # 다운로드 버튼 (오른쪽 정렬)
            st.markdown("---")
            col_left, col_right = st.columns([3, 1])
            with col_right:
                csv_calc = results['calc_df'].to_csv(index=False)
                st.download_button(
                    label="📥 계산된 곡선 다운로드",
                    data=csv_calc,
                    file_name="MM_calculated_curves.csv",
                    mime="text/csv",
                    use_container_width=True
                )
        
        with tab2:
            st.subheader("Michaelis-Menten Fitting Results")
            # Concentration과 Concentration [ug/mL] 열 제외하고 표시, 열 이름에 단위 추가
            display_df = results['results_df'].drop(columns=['Concentration', 'Concentration [ug/mL]'], errors='ignore').copy()
            display_df = display_df.rename(columns={
                'Vmax': 'Vmax [RFU]',
                'Km': 'Km [min]',
                'F0': 'F0 [RFU]',
                'Fmax': 'Fmax [RFU]',
                'k': 'k [min⁻¹]',
                'R_squared': 'R²'
            })
            # 열 순서 지정: Type, Equation, Vmax, Km, F0, Fmax, k, R²
            column_order = ['Type', 'Equation', 'Vmax [RFU]', 'Km [min]', 'F0 [RFU]', 'Fmax [RFU]', 'k [min⁻¹]', 'R²']
            display_df = display_df[column_order]
            st.dataframe(display_df, use_container_width=True, hide_index=True)
            
            # 다운로드 버튼 (오른쪽 정렬)
            st.markdown("---")
            col_left, col_right = st.columns([3, 1])
            with col_right:
                csv_results = results['results_df'].to_csv(index=False)
                st.download_button(
                    label="📥 결과 다운로드",
                    data=csv_results,
                    file_name="MM_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )
        
        with tab3:
            st.subheader("추가 데이터 다운로드")
            st.caption("Calibration Curve 및 Fit Curves 원본 데이터")
            
            # CSV 다운로드 버튼들
            col1, col2 = st.columns(2)
            
            with col1:
                csv_cal = results['cal_curve_df'].to_csv(index=False)
                st.download_button(
                    label="📥 Calibration Curve (CSV)",
                    data=csv_cal,
                    file_name="MM_calibration_curve.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            with col2:
                csv_fits = results['fit_curves_df'].to_csv(index=False)
                st.download_button(
                    label="📥 Fit Curves (CSV)",
                    data=csv_fits,
                    file_name="MM_fit_curves.csv",
                    mime="text/csv",
                    use_container_width=True
                )



