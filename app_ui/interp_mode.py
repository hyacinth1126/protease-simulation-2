import os
import numpy as np
import pandas as pd
import streamlit as st

from data_interpolation_mode.interpolate_prism import (
    exponential_association,
    create_prism_interpolation_range
)


def data_interpolation_mode(st):
    """Data Interpolation 모드 - Prism 스타일 보간"""
    
    # 폴더 구조 생성
    os.makedirs("prep_data/raw", exist_ok=True)
    os.makedirs("data_interpolation_mode/results", exist_ok=True)
    
    st.header("📈 Data Interpolation 모드")
    st.markdown("GraphPad스타일 보간 - Fitting 결과에서 고밀도 곡선 생성")
    st.markdown("---")
    
    # 사이드바 설정
    st.sidebar.title("⚙️ Data Interpolation 설정")
    
    # 데이터 업로드
    st.sidebar.subheader("📁 데이터 업로드")
    
    # MM Results 파일 업로드
    mm_file = st.sidebar.file_uploader(
        "MM Results CSV 파일 업로드",
        type=['csv'],
        help="MM_results_detailed.csv: Fitting 파라미터 포함",
        key="mm_results_upload"
    )
    
    # Raw data 파일 업로드 (시간 범위 확인용)
    raw_file = st.sidebar.file_uploader(
        "Prep Raw CSV 파일 업로드 (선택사항)",
        type=['csv'],
        help="prep_raw.csv: 시간 범위 확인용",
        key="raw_data_upload"
    )
    
    # 샘플 데이터 다운로드
    col1, col2 = st.sidebar.columns(2)
    with col1:
        try:
            with open("prep_raw_data_mode/results/MM_results_detailed.csv", "rb") as f:
                sample_bytes = f.read()
            st.download_button(
                label="📥 샘플 MM Results",
                data=sample_bytes,
                file_name="MM_results_sample.csv",
                mime="text/csv"
            )
        except Exception:
            pass
    
    with col2:
        try:
            with open("prep_data/raw/prep_raw.csv", "rb") as f:
                sample_bytes = f.read()
            st.download_button(
                label="📥 샘플 Raw Data",
                data=sample_bytes,
                file_name="prep_raw_sample.csv",
                mime="text/csv"
            )
        except Exception:
            pass
    
    # 보간 설정
    st.sidebar.markdown("---")
    st.sidebar.subheader("⚙️ 보간 설정")
    
    n_points = st.sidebar.slider(
        "보간 포인트 개수",
        min_value=100,
        max_value=2000,
        value=1000,
        step=100,
        help="더 많은 포인트 = 더 부드러운 곡선"
    )
    
    include_y_to_x = st.sidebar.checkbox(
        "Y → X 역보간 포함",
        value=False,
        help="특정 RFU 값에 대한 시간 계산"
    )
    
    # 데이터 로드
    mm_results_df = None
    raw_data_df = None
    
    # MM Results 읽기
    if mm_file is not None:
        import tempfile
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv', mode='wb') as tmp_file:
            tmp_file.write(mm_file.getbuffer())
            tmp_path = tmp_file.name
        
        try:
            mm_results_df = pd.read_csv(tmp_path)
            os.unlink(tmp_path)
        except Exception as e:
            st.error(f"MM Results 파일 읽기 오류: {e}")
            os.unlink(tmp_path)
            return
    else:
        # 기본 샘플 데이터 사용
        try:
            mm_results_df = pd.read_csv('prep_raw_data_mode/results/MM_results_detailed.csv')
            st.sidebar.info("prep_raw_data_mode/results/MM_results_detailed.csv 사용 중")
        except FileNotFoundError:
            st.error("MM Results 파일을 찾을 수 없습니다. CSV 파일을 업로드해주세요.")
            st.info("💡 먼저 'Prep Raw Data 모드'에서 분석을 실행하여 MM_results_detailed.csv를 생성하거나, 파일을 업로드해주세요.")
            return
    
    # Raw data 읽기 (선택사항)
    if raw_file is not None:
        import tempfile
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv', mode='wb') as tmp_file:
            tmp_file.write(raw_file.getbuffer())
            tmp_path = tmp_file.name
        
        try:
            raw_data_df = pd.read_csv(tmp_path, sep='\t', skiprows=[0, 1])
            os.unlink(tmp_path)
        except Exception as e:
            st.warning(f"Raw data 파일 읽기 오류 (무시됨): {e}")
            os.unlink(tmp_path)
    else:
        try:
            raw_data_df = pd.read_csv('prep_data/raw/prep_raw.csv', sep='\t', skiprows=[0, 1])
        except Exception:
            pass
    
    # 데이터 미리보기
    st.subheader("📋 MM Fitting Results 미리보기")
    
    if 'Concentration [ug/mL]' in mm_results_df.columns:
        st.metric("농도 조건 수", len(mm_results_df))
        
        with st.expander("📊 MM Fitting Parameters"):
            display_cols = ['Concentration', 'Concentration [ug/mL]', 'F0', 'Fmax', 'k', 'R_squared']
            available_cols = [col for col in display_cols if col in mm_results_df.columns]
            st.dataframe(mm_results_df[available_cols], use_container_width=True, height=300)
    else:
        st.dataframe(mm_results_df, use_container_width=True, height=300)
    
    # 시간 범위 확인
    x_data_min = 0
    x_data_max = 30  # 기본값
    
    if raw_data_df is not None:
        try:
            time_col = raw_data_df.columns[0]
            times = pd.to_numeric(raw_data_df[time_col].values, errors='coerce')
            times = times[~np.isnan(times)]
            if len(times) > 0:
                x_data_min = float(np.min(times))
                x_data_max = float(np.max(times))
                st.info(f"📊 데이터 시간 범위: {x_data_min:.1f} - {x_data_max:.1f} min")
        except Exception:
            st.warning("Raw data에서 시간 범위를 추출할 수 없습니다. 기본값(0-30 min) 사용")
    
    # 보간 실행
    st.markdown("---")
    
    if st.button("🚀 Prism 스타일 보간 실행", type="primary", use_container_width=True):
        with st.spinner("보간 중..."):
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # 1. 보간 범위 계산
            status_text.text("1️⃣ 보간 범위 계산 중...")
            
            times_array = np.array([x_data_min, x_data_max])
            x_range_min, x_range_max = create_prism_interpolation_range(times_array)
            
            st.info(f"📐 보간 범위: {x_range_min:.3f} - {x_range_max:.3f} min (데이터: {x_data_min:.1f} - {x_data_max:.1f} min)")
            
            progress_bar.progress(0.2)
            
            # 2. X → Y 보간 수행
            status_text.text("2️⃣ X → Y 보간 수행 중...")
            
            x_interp = np.linspace(x_range_min, x_range_max, n_points + 1)
            
            all_interp_data = []
            
            for idx, row in mm_results_df.iterrows():
                conc_name = row.get('Concentration', f'Conc_{idx}')
                F0 = row['F0']
                Fmax = row['Fmax']
                k = row['k']
                
                # X → Y 보간
                y_interp = exponential_association(x_interp, F0, Fmax, k)
                
                for x, y in zip(x_interp, y_interp):
                    all_interp_data.append({
                        'Concentration': conc_name,
                        'Concentration [ug/mL]': row.get('Concentration [ug/mL]', None),
                        'Time_min': x,
                        'RFU_Interpolated': y,
                        'Is_Extrapolated': (x < x_data_min) or (x > x_data_max)
                    })
            
            interp_df = pd.DataFrame(all_interp_data)
            
            progress_bar.progress(0.6)
            
            # 3. Y → X 역보간 (선택사항)
            y_to_x_df = None
            
            if include_y_to_x:
                status_text.text("3️⃣ Y → X 역보간 수행 중...")
                
                y_to_x_examples = []
                
                for idx, row in mm_results_df.iterrows():
                    conc_name = row.get('Concentration', f'Conc_{idx}')
                    F0 = row['F0']
                    Fmax = row['Fmax']
                    k = row['k']
                    
                    # Y 값 예제 (F0에서 Fmax까지 5개)
                    if Fmax > F0:
                        y_examples = np.linspace(F0 + (Fmax - F0) * 0.1, 
                                                Fmax - (Fmax - F0) * 0.1, 5)
                        
                        for y in y_examples:
                            # 역함수로 X 계산
                            try:
                                if k > 0:
                                    x_calc = -np.log(1 - (y - F0) / (Fmax - F0)) / k
                                    
                                    if x_range_min <= x_calc <= x_range_max:
                                        y_to_x_examples.append({
                                            'Concentration': conc_name,
                                            'Target_RFU': y,
                                            'Calculated_Time_min': x_calc,
                                            'Is_In_Data_Range': (x_data_min <= x_calc <= x_data_max)
                                        })
                            except Exception:
                                continue
                
                if y_to_x_examples:
                    y_to_x_df = pd.DataFrame(y_to_x_examples)
            
            progress_bar.progress(0.9)
            
            # 4. 결과 저장
            status_text.text("4️⃣ 결과 저장 중...")
            
            # 결과 파일을 interpolation_results 폴더에 자동 저장
            try:
                # Interpolated curves 저장
                interp_df.to_csv('data_interpolation_mode/results/MM_interpolated_curves.csv', index=False)
                
                # Y to X results 저장 (있을 경우)
                if y_to_x_df is not None:
                    y_to_x_df.to_csv('data_interpolation_mode/results/MM_Y_to_X_interpolation.csv', index=False)
                
                st.sidebar.success("✅ 결과 파일이 data_interpolation_mode/results/ 에 저장되었습니다!")
            except Exception as e:
                st.sidebar.warning(f"⚠️ 파일 저장 중 오류: {e}")
            
            st.session_state['interpolation_results'] = {
                'interp_df': interp_df,
                'y_to_x_df': y_to_x_df,
                'x_range_min': x_range_min,
                'x_range_max': x_range_max,
                'x_data_min': x_data_min,
                'x_data_max': x_data_max,
                'mm_results_df': mm_results_df
            }
            
            progress_bar.progress(1.0)
            status_text.text("✅ 보간 완료!")
    
    # 결과 표시
    if 'interpolation_results' in st.session_state:
        results = st.session_state['interpolation_results']
        
        st.markdown("---")
        st.subheader("📊 보간 결과")
        
        # 탭 구성
        tabs = ["📈 Interpolated Curves", "📋 Data Table", "💾 Download"]
        if results['y_to_x_df'] is not None:
            tabs.insert(2, "🔄 Y → X Results")
        
        tab_objects = st.tabs(tabs)
        
        # Tab 1: 그래프
        with tab_objects[0]:
            st.subheader("Interpolated Curves")
            
            import plotly.graph_objects as go
            fig = go.Figure()
            
            colors = ['blue', 'red', 'orange', 'green', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
            
            # 농도 순서대로 정렬
            if 'Concentration [ug/mL]' in results['mm_results_df'].columns:
                conc_order = results['mm_results_df'].sort_values('Concentration [ug/mL]')['Concentration'].tolist()
            else:
                conc_order = results['mm_results_df']['Concentration'].tolist()
            
            x_data_min = results['x_data_min']
            x_data_max = results['x_data_max']
            
            for idx, conc_name in enumerate(conc_order):
                color = colors[idx % len(colors)]
                
                # 보간 곡선
                subset = results['interp_df'][results['interp_df']['Concentration'] == conc_name]
                
                # 데이터 범위 내 (실선)
                interp_in_range = subset[~subset['Is_Extrapolated']]
                if len(interp_in_range) > 0:
                    fig.add_trace(go.Scatter(
                        x=interp_in_range['Time_min'],
                        y=interp_in_range['RFU_Interpolated'],
                        mode='lines',
                        name=f'{conc_name} (Interpolated)',
                        line=dict(color=color, width=2.5),
                        legendgroup=conc_name,
                        showlegend=True
                    ))
                
                # 외삽 영역 (점선)
                interp_extrap = subset[subset['Is_Extrapolated']]
                if len(interp_extrap) > 0:
                    fig.add_trace(go.Scatter(
                        x=interp_extrap['Time_min'],
                        y=interp_extrap['RFU_Interpolated'],
                        mode='lines',
                        name=f'{conc_name} (Extrapolated)',
                        line=dict(color=color, width=2, dash='dash'),
                        opacity=0.5,
                        legendgroup=conc_name,
                        showlegend=False
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
            
            fig.update_xaxes(range=[results['x_range_min'] - 1, results['x_range_max'] + 1])
            fig.update_yaxes(rangemode='tozero')
            
            st.plotly_chart(fig, use_container_width=True)
        
        # Tab 2: 데이터 테이블
        with tab_objects[1]:
            st.subheader("Interpolation Data")
            
            # 요약 통계
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("총 포인트 수", len(results['interp_df']))
            with col2:
                in_range = results['interp_df'][~results['interp_df']['Is_Extrapolated']]
                st.metric("보간 포인트", len(in_range))
            with col3:
                extrap = results['interp_df'][results['interp_df']['Is_Extrapolated']]
                st.metric("외삽 포인트", len(extrap))
            
            st.markdown("---")
            
            # 데이터 테이블
            st.dataframe(results['interp_df'], use_container_width=True, height=400)
        
        # Tab 3: Y → X 결과 (선택사항)
        if results['y_to_x_df'] is not None:
            with tab_objects[2]:
                st.subheader("Y → X Interpolation Results")
                st.caption("특정 RFU 값에 도달하는 시간 계산")
                
                st.dataframe(results['y_to_x_df'], use_container_width=True, height=400)
        
        # Download 탭
        download_tab_idx = 3 if results['y_to_x_df'] is not None else 2
        with tab_objects[download_tab_idx]:
            st.subheader("데이터 다운로드")
            
            col1, col2 = st.columns(2)
            
            with col1:
                csv_interp = results['interp_df'].to_csv(index=False)
                st.download_button(
                    label="📥 Interpolated Data (CSV)",
                    data=csv_interp,
                    file_name="MM_interpolated_curves.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            if results['y_to_x_df'] is not None:
                with col2:
                    csv_y_to_x = results['y_to_x_df'].to_csv(index=False)
                    st.download_button(
                        label="📥 Y → X Results (CSV)",
                        data=csv_y_to_x,
                        file_name="MM_Y_to_X_interpolation.csv",
                        mime="text/csv",
                        use_container_width=True
                    )



