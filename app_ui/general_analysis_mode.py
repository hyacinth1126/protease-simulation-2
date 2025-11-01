import pandas as pd
import streamlit as st

from general_analysis_mode.analysis import (
    UnitStandardizer,
    DataNormalizer,
    RegionDivider,
    ModelA_SubstrateDepletion,
    ModelB_EnzymeDeactivation,
    ModelC_MassTransfer,
    ModelD_ConcentrationDependentFmax,
    ModelE_ProductInhibition,
    ModelF_EnzymeSurfaceSequestration
)
from general_analysis_mode.plot import Visualizer


def verbose_callback(message: str, level: str = "info"):
    """Callback function for logging from analysis modules"""
    if level == "error":
        st.error(message)
    elif level == "warning":
        st.warning(message)
    elif level == "debug":
        st.code(message)
    else:
        st.info(message)


def general_analysis_mode(st):
    """일반 분석 모드 - 표준 FRET 분석"""
    
    # Sidebar configuration
    st.sidebar.title("⚙️ 설정")
    
    enzyme_mw = st.sidebar.number_input(
        "효소 분자량 (kDa)",
        min_value=1.0,
        max_value=500.0,
        value=56.6,
        step=0.1,
        help="Kgp: 56.6 kDa"
    )
    
    enzyme_name = st.sidebar.text_input(
        "효소 이름 (선택사항)",
        value="",
        placeholder="enzyme",
        help="그래프 범례에 표시될 효소 이름 (비워두면 'enzyme' 표시)"
    )
    if enzyme_name.strip() == "":
        enzyme_name = "enzyme"
    
    substrate_name = st.sidebar.text_input(
        "기질 이름 (선택사항)",
        value="",
        placeholder="substrate",
        help="그래프 범례에 표시될 기질 이름 (비워두면 'substrate' 표시)"
    )
    if substrate_name.strip() == "":
        substrate_name = "substrate"
    # 구분선 후 데이터 소스 섹션
    st.sidebar.markdown("---")
    st.sidebar.subheader("📁 데이터 소스")
    
    # 데이터 소스 타입 선택
    data_source_type = st.sidebar.radio(
        "데이터 소스 선택",
        ["Raw Data Points", "Fitted Curves (from Prep mode)"],
        help="Raw Data Points: 원본 측정 데이터 | Fitted Curves: Prep 모드에서 생성된 fitting/interpolation 곡선"
    )

    uploaded_file = st.sidebar.file_uploader(
        "CSV 파일 업로드",
        type=['csv'],
        help="Raw Data: time_s, enzyme_ugml, FL_intensity, SD | Fitted Curves: Concentration, Time_min, RFU_*"
    )
    # Provide sample data download based on data source type
    if data_source_type == "Raw Data Points":
        try:
            with open("fitc_peptide_timeseries.csv", "rb") as f:
                sample_bytes = f.read()
            st.sidebar.download_button(
                label="샘플 원본 데이터 다운로드 (CSV)",
                data=sample_bytes,
                file_name="raw_data.csv",
                mime="text/csv",
                help="배포된 기본 CSV를 다운로드합니다."
            )
        except Exception:
            pass
    else:
        # Fitted Curves 샘플 다운로드
        col1, col2 = st.sidebar.columns(2)
        with col1:
            try:
                with open("prep_raw_data_mode/results/MM_calculated_curves.csv", "rb") as f:
                    sample_bytes = f.read()
                st.download_button(
                    label="📥 Calculated",
                    data=sample_bytes,
                    file_name="calculated_curves_sample.csv",
                    mime="text/csv"
                )
            except Exception:
                pass
        with col2:
            try:
                with open("data_interpolation_mode/results/MM_interpolated_curves.csv", "rb") as f:
                    sample_bytes = f.read()
                st.download_button(
                    label="📥 Interpolated",
                    data=sample_bytes,
                    file_name="interpolated_curves_sample.csv",
                    mime="text/csv"
                )
            except Exception:
                pass
    
    # Step 1: Load data based on source type
    if data_source_type == "Raw Data Points":
        # Raw Data Points 모드
        if uploaded_file is not None:
            df_raw = pd.read_csv(uploaded_file)
        else:
            # Use default sample data
            try:
                df_raw = pd.read_csv("fitc_peptide_timeseries.csv")
                st.sidebar.info("fitc_peptide_timeseries.csv 사용 중")
            except FileNotFoundError:
                st.error("데이터 파일을 찾을 수 없습니다. CSV 파일을 업로드해주세요.")
                st.stop()
    else:
        # Fitted Curves 모드
        if uploaded_file is not None:
            df_fitted = pd.read_csv(uploaded_file)
        else:
            # Try to load from fitting_results first, then interpolation_results
            try:
                df_fitted = pd.read_csv("prep_raw_data_mode/results/MM_calculated_curves.csv")
                st.sidebar.info("prep_raw_data_mode/results/MM_calculated_curves.csv 사용 중")
            except FileNotFoundError:
                try:
                    df_fitted = pd.read_csv("data_interpolation_mode/results/MM_interpolated_curves.csv")
                    st.sidebar.info("data_interpolation_mode/results/MM_interpolated_curves.csv 사용 중")
                except FileNotFoundError:
                    st.error("Fitted curves 파일을 찾을 수 없습니다. 먼저 'Prep Raw Data 모드' 또는 'Data Interpolation 모드'를 실행하거나 파일을 업로드해주세요.")
                    st.stop()
        
        # Convert fitted curves to raw data format
        # Fitted curves format: Concentration, Concentration [ug/mL], Time_min, RFU_Calculated/RFU_Interpolated, Is_Extrapolated
        # Target format: time_min, enzyme_ugml, FL_intensity, SD
        
        # Detect RFU column name
        rfu_col = None
        if 'RFU_Calculated' in df_fitted.columns:
            rfu_col = 'RFU_Calculated'
        elif 'RFU_Interpolated' in df_fitted.columns:
            rfu_col = 'RFU_Interpolated'
        else:
            st.error("RFU 데이터 컬럼을 찾을 수 없습니다. (RFU_Calculated 또는 RFU_Interpolated)")
            st.stop()
        
        # Pivot data to get one row per time point with columns for each concentration
        df_raw_converted = []
        
        # Get unique time points
        unique_times = sorted(df_fitted['Time_min'].unique())
        
        for time in unique_times:
            time_data = df_fitted[df_fitted['Time_min'] == time]
            
            # Create row for each concentration
            for _, row in time_data.iterrows():
                conc_ugml = row.get('Concentration [ug/mL]', 0)
                rfu = row[rfu_col]
                
                df_raw_converted.append({
                    'time_min': time,
                    'enzyme_ugml': conc_ugml,
                    'FL_intensity': rfu,
                    'SD': 0  # Fitted curves don't have SD
                })
        
        df_raw = pd.DataFrame(df_raw_converted)
        
        st.sidebar.success(f"✅ {len(df_fitted['Concentration'].unique())}개 농도 조건, {len(unique_times)}개 시간 포인트 로드됨")
        st.sidebar.info(f"📊 RFU 컬럼: {rfu_col}")
    
    # Store data source type for later use
    st.session_state['data_source_type'] = data_source_type
    
    # Step 2: Standardize units
    standardizer = UnitStandardizer(enzyme_mw=enzyme_mw)
    df_standardized = standardizer.standardize(df_raw)
    
    # Store time unit for later use
    time_unit = 'min' if 'time_min' in df_raw.columns else 's'
    st.session_state['time_unit'] = time_unit
    
    # Step 3-4: Iterative normalization and region division
    normalizer = DataNormalizer()
    region_divider = RegionDivider()
    
    # Read iteration setting from session (set in 정규화 탭); default 2
    max_iterations = int(st.session_state.get('max_iterations', 2))
    
    # Step 3-1: Initial temporary normalization (model-free threshold)
    df_current = normalizer.normalize_temporary(df_standardized)
    
    # Iterative loop: Divide regions → Final normalization → Divide regions → ...
    for iteration in range(max_iterations):
        
        # Step 4: Divide regions
        df_current = region_divider.divide_regions(df_current)
        
        # Step 3-2: Final normalization (using current region information)
        df_current = normalizer.normalize_final(df_current)
    
    df = df_current
    
    # Display data
    st.subheader("📊 데이터 미리보기")
    
    # Detect original column names for display
    time_unit = st.session_state.get('time_unit', 's')
    if time_unit == 'min':
        time_display = f"0 - {df['time_s'].max():.0f} 분"
        time_label = "시간 (분)"
    else:
        time_display = f"0 - {df['time_s'].max():.0f} 초" if df['time_s'].max() < 100 else f"0 - {df['time_s'].max()/60:.1f} 분"
        time_label = "시간 (초)"
    # Determine concentration unit from normalized data
    conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
    if 'uM' in conc_col:
        conc_unit = "μM"
    elif 'nM' in conc_col:
        conc_unit = "nM"
    else:
        conc_unit = "μg/mL"
    
    st.session_state['time_label'] = time_label
    st.session_state['conc_unit'] = conc_unit
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("데이터 포인트", len(df))
    with col2:
        st.metric(f"농도 조건 ({conc_unit})", df[conc_col].nunique())
    with col3:
        st.metric("시간 범위", time_display)
    
    # Tabs for different views
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📈 원본 데이터", 
        "📊 정규화 데이터", 
        "🔬 모델 피팅",
        "📉 모델 비교",
        "💡 진단 분석"
    ])
    
    with tab1:
        st.plotly_chart(
            Visualizer.plot_raw_data(df, conc_unit, time_label, 
                                    use_lines=(st.session_state.get('data_source_type') == 'Fitted Curves (from Prep mode)'),
                                    enzyme_name=enzyme_name, 
                                    substrate_name=substrate_name), 
            use_container_width=True
        )
        
        st.subheader("Raw data table")
        st.dataframe(df, height=400, use_container_width=True)
    
    with tab2:
        # Controls and method description for normalization
        st.subheader("정규화 설정 및 방법")
        st.caption("최종 정규화와 구간 구분을 반복하며 수렴시킵니다.")
        st.number_input(
            "정규화-구간 반복 횟수",
            min_value=2,
            max_value=10,
            value=int(st.session_state.get('max_iterations', 2)),
            step=1,
            key="max_iterations",
            help="최소 2회 이상 권장. 값을 변경하면 화면이 다시 계산됩니다."
        )
        with st.expander("정규화 방법 보기", expanded=False):
            st.markdown("""
            - 각 농도별 지수 피팅: F(t) = F₀ + A·(1−e⁻ᵏᵗ)
            - 점근선 Fmax = F₀ + A 사용
            - α(t) = (F(t) − F₀) / (Fmax − F₀)
            """)
        st.markdown(f"현재 반복 횟수: **{int(st.session_state.get('max_iterations', 2))}**")

        st.plotly_chart(
            Visualizer.plot_normalized_data(df, conc_unit, time_label, 
                                           use_lines=(st.session_state.get('data_source_type') == 'Fitted Curves (from Prep mode)'),
                                           enzyme_name=enzyme_name,
                                           substrate_name=substrate_name), 
            use_container_width=True
        )
        
        # Summary statistics
        st.subheader("정규화 요약 (지수 피팅 기반)")
        
        summary_data = []
        for conc in sorted(df[conc_col].unique()):
            subset = df[df[conc_col] == conc]
            # Check if optional columns exist
            fmax_std = f"{subset['Fmax_std'].iloc[0]:.1f}" if 'Fmax_std' in subset.columns else "N/A"
            fit_k = f"{subset['fit_k'].iloc[0]:.4f}" if 'fit_k' in subset.columns else "N/A"
            
            summary_data.append({
                f'농도 ({conc_unit})': conc,
                'F0 (초기)': f"{subset['F0'].iloc[0]:.1f}",
                'Fmax (점근선)': f"{subset['Fmax'].iloc[0]:.1f}",
                'Fmax 표준편차': fmax_std,
                '피팅 k (s⁻¹)': fit_k,
                'α 범위': f"{subset['alpha'].min():.3f} - {subset['alpha'].max():.3f}",
                'α 평균': f"{subset['alpha'].mean():.3f}"
            })
        
        summary_df = pd.DataFrame(summary_data)
        st.dataframe(summary_df, use_container_width=True)
        
        st.info("📊 각 농도별로 F(t) = F0 + A·(1-exp(-k·t)) 형태의 지수 함수를 피팅하여 점근선 Fmax를 결정합니다.")
    
    with tab3:
        st.subheader("🔬 글로벌 모델 피팅")
        
        st.markdown("""
        **기본 모델 (A-C)**: 고전적 효소 키네틱 메커니즘  
        **확장 모델 (D-F)**: Fmax 농도 의존성 설명 (겔 침투, 생성물 억제, 효소 흡착)
        """)
        
        # Model selection
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**기본 모델**")
            fit_model_a = st.checkbox("모델 A: 기질 고갈", value=True)
            st.caption("✓ 1차 반응 및 기질 고갈")
            
            fit_model_b = st.checkbox("모델 B: 효소 비활성화", value=True)
            st.caption("✓ 효소 비활성화 & 시간 의존")
            
            fit_model_c = st.checkbox("모델 C: 물질전달 제한", value=True)
            st.caption("✓ 확산 제한 & 접근성 제약")
        
        with col2:
            st.markdown("**확장 모델 (Fmax 의존성)**")
            fit_model_d = st.checkbox("모델 D: 농도 의존 Fmax", value=True)
            st.caption("✓ 겔 침투 깊이 & 2차 절단")
            
            fit_model_e = st.checkbox("모델 E: 생성물 억제", value=True)
            st.caption("✓ 생성물 축적 & 경쟁 억제")
            
            fit_model_f = st.checkbox("모델 F: 효소 흡착/격리", value=True)
            st.caption("✓ 표면 흡착 & 비가역 결합")
        
        if st.button("🚀 글로벌 피팅 실행", type="primary"):
            results = []
            
            # Create a status container
            status_container = st.empty()
            result_container = st.container()
            
            # Model A
            if fit_model_a:
                with status_container:
                    with st.spinner("🔄 모델 A: 기질 고갈 피팅 중..."):
                        model_a = ModelA_SubstrateDepletion(enzyme_mw=enzyme_mw)
                        result_a = model_a.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_a)
                
                if result_a:
                    with result_container:
                        st.success(f"✅ 모델 A 완료: R² = {result_a.r_squared:.4f}, AIC = {result_a.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ 모델 A 피팅 실패")
            
            # Model B
            if fit_model_b:
                with status_container:
                    with st.spinner("🔄 모델 B: 효소 비활성화 피팅 중..."):
                        model_b = ModelB_EnzymeDeactivation(enzyme_mw=enzyme_mw)
                        result_b = model_b.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_b)
                
                if result_b:
                    with result_container:
                        st.success(f"✅ 모델 B 완료: R² = {result_b.r_squared:.4f}, AIC = {result_b.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ 모델 B 피팅 실패")
            
            # Model C
            if fit_model_c:
                with status_container:
                    with st.spinner("🔄 모델 C: 물질전달 제한 피팅 중..."):
                        model_c = ModelC_MassTransfer(enzyme_mw=enzyme_mw)
                        result_c = model_c.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_c)
                
                if result_c:
                    with result_container:
                        st.success(f"✅ 모델 C 완료: R² = {result_c.r_squared:.4f}, AIC = {result_c.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ 모델 C 피팅 실패")
            
            # Model D
            if fit_model_d:
                with status_container:
                    with st.spinner("🔄 모델 D: 농도 의존 Fmax 피팅 중..."):
                        model_d = ModelD_ConcentrationDependentFmax(enzyme_mw=enzyme_mw)
                        result_d = model_d.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_d)
                
                if result_d:
                    with result_container:
                        st.success(f"✅ 모델 D 완료: R² = {result_d.r_squared:.4f}, AIC = {result_d.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ 모델 D 피팅 실패")
            
            # Model E
            if fit_model_e:
                with status_container:
                    with st.spinner("🔄 모델 E: 생성물 억제 피팅 중..."):
                        model_e = ModelE_ProductInhibition(enzyme_mw=enzyme_mw)
                        result_e = model_e.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_e)
                
                if result_e:
                    with result_container:
                        st.success(f"✅ 모델 E 완료: R² = {result_e.r_squared:.4f}, AIC = {result_e.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ 모델 E 피팅 실패")
            
            # Model F
            if fit_model_f:
                with status_container:
                    with st.spinner("🔄 모델 F: 효소 흡착/격리 피팅 중..."):
                        model_f = ModelF_EnzymeSurfaceSequestration(enzyme_mw=enzyme_mw)
                        result_f = model_f.fit_global(df, verbose_callback=verbose_callback)
                        results.append(result_f)
                
                if result_f:
                    with result_container:
                        st.success(f"✅ 모델 F 완료: R² = {result_f.r_squared:.4f}, AIC = {result_f.aic:.2f}")
                else:
                    with result_container:
                        st.error("❌ 모델 F 피팅 실패")
            
            # Clear status container after all done
            status_container.empty()
            
            # Store results in session state
            st.session_state['fit_results'] = results
            st.session_state['df'] = df
            
            # Show completion message
            with result_container:
                st.success("🎉 모든 모델 피팅 완료! '모델 비교' 탭에서 결과를 확인하세요.")
    
    with tab4:
        if 'fit_results' in st.session_state:
            results = st.session_state['fit_results']
            df = st.session_state['df']
            
            st.subheader("📊 모델 비교")
            
            # Comparison table
            comparison_df = Visualizer.create_comparison_table(results)
            st.dataframe(comparison_df, use_container_width=True)
            
            # Determine best model
            valid_results = [r for r in results if r is not None]
            if valid_results:
                best_aic = min(r.aic for r in valid_results)
                best_model = [r for r in valid_results if r.aic == best_aic][0]
                
                st.success(f"🏆 최적 모델 (최저 AIC): **{best_model.name}** (AIC = {best_model.aic:.2f})")
                
                # Parameter details for best model
                st.subheader(f"최적 모델 파라미터: {best_model.name}")
                param_data = []
                for param, value in best_model.params.items():
                    std = best_model.params_std.get(param, 0)
                    param_data.append({
                        '파라미터': param,
                        '값': f"{value:.4e}",
                        '표준오차': f"{std:.4e}",
                        '상대오차': f"{(std/value*100):.2f}%" if value != 0 else "N/A"
                    })
                st.dataframe(pd.DataFrame(param_data), use_container_width=True)
            
            # Plot all model fits
            st.subheader("📈 전체 모델 피팅 결과")
            st.plotly_chart(
                Visualizer.plot_model_fits(df, results, conc_unit, time_label,
                                          enzyme_name=enzyme_name,
                                          substrate_name=substrate_name), 
                use_container_width=True
            )
            
            # Individual model plots
            st.subheader("📊 개별 모델 비교")
            st.markdown("각 모델별로 원본 데이터와 피팅 결과를 비교합니다.")
            
            # Create tabs for each model
            model_names = [r.name for r in results if r is not None]
            
            if len(model_names) > 0:
                model_tabs_ui = st.tabs(model_names)
                
                for idx, (tab, result) in enumerate(zip(model_tabs_ui, [r for r in results if r is not None])):
                    with tab:
                        # Color scheme for each model
                        model_colors = ['#FF6B6B', '#4ECDC4', '#FFD93D']
                        color = model_colors[idx % len(model_colors)]
                        
                        # Display individual model plot
                        st.plotly_chart(
                            Visualizer.plot_individual_model(df, result, conc_unit, time_label, color),
                            use_container_width=True
                        )
                        
                        # Display parameters
                        st.markdown(f"**{result.name} 파라미터**")
                        param_cols = st.columns(len(result.params))
                        for col_idx, (param, value) in enumerate(result.params.items()):
                            with param_cols[col_idx]:
                                std = result.params_std.get(param, 0)
                                st.metric(
                                    label=param,
                                    value=f"{value:.4e}",
                                    delta=f"±{std:.4e}" if std > 0 else None
                                )
            
            # Download results
            st.subheader("💾 결과 다운로드")
            csv = comparison_df.to_csv(index=False)
            st.download_button(
                label="비교 테이블 다운로드 (CSV)",
                data=csv,
                file_name="model_comparison.csv",
                mime="text/csv"
            )
        else:
            st.info("👈 먼저 '모델 피팅' 탭에서 피팅을 실행해주세요.")
    
    with tab5:
        st.subheader("💡 진단 분석")
        
        # Initial rate analysis
        st.plotly_chart(
            Visualizer.plot_initial_rates(df, conc_unit, time_unit), 
            use_container_width=True
        )
        
        st.markdown("""
        ### 📋 모델 선택 가이드라인
        
        #### 기본 모델 (A-C)
        
        **모델 A (기질 고갈)** 선호 조건:
        - 초기 속도 v₀가 [E]에 대해 선형 관계 (낮은 [E]에서)
        - 포화 형광 F∞ ≈ 일정 (정규화된 α → 1)
        - 유의미한 효소 비활성화가 관찰되지 않음
        
        **모델 B (효소 비활성화)** 선호 조건:
        - F∞ < 이론적 최대값 (포화에서 α < 1)
        - 빠른 초기 증가 후 예상보다 낮은 수준에서 평탄화
        - kd > 0이며 유의미한 기여도
        
        **모델 C (물질전달 제한)** 선호 조건:
        - 초기 버스트(0-5초) 후 느린 접근
        - 교반/유속에 민감
        - 높은 [E]에서 v₀ vs [E] 그래프가 포화 양상
        
        #### 확장 모델 (D-F): **Fmax가 [E]에 따라 변하는 경우**
        
        **모델 D (농도 의존 Fmax)** 선호 조건:
        - 높은 [E]에서 α_max 증가 (더 많은 기질 접근)
        - 겔 침투 깊이 효과 (두꺼운/밀집 겔)
        - 2차 절단으로 생성물 방출 증가
        - **파라미터**: α_∞ (최대값), k_access (접근성 계수)
        
        **모델 E (생성물 억제)** 선호 조건:
        - 초기 빠른 증가 후 감속 (생성물 축적)
        - 낮은 [E]에서 더 큰 억제 효과
        - 생성물 제거 시 반응 속도 회복
        - **파라미터**: Ki_eff (억제 상수)
        
        **모델 F (효소 흡착/격리)** 선호 조건:
        - 높은 [E]에서 상대적으로 덜 영향받음 (포화)
        - 음전하 표면/PDA 코팅, 밀집 겔 구조
        - 시간에 따른 효소 활성 감소 (비가역)
        - **파라미터**: k_ads (흡착속도), K_ads (평형상수)
        
        ### 📊 통계 기준
        - **AIC/BIC**: 낮을수록 좋음 (파라미터 수 페널티)
        - **R²**: 높을수록 좋음 (>0.95 우수)
        - **RMSE**: 낮을수록 좋음
        - **Δ AIC > 10**: 높은 AIC 모델에 대한 강력한 반증
        - **Δ AIC < 2**: 모델 간 유의미한 차이 없음
        """)
        
        # Experimental suggestions
        st.subheader("🧪 제안 후속 실험 (모델 구분)")
        
        st.markdown("""
        ### 🔍 Fmax 농도 의존성 확인 실험
        
        1. **다양한 [E]에서 장시간 측정** (30분-1시간)
           - 각 농도별 포화 형광값(Fmax) 정량 측정
           - [E] vs Fmax 플롯 → 선형/포화 양상 확인
           - **선형 증가** → 모델 D 가능성
           - **일정** → 기본 모델 A-C
        
        2. **겔 두께 변화 테스트** (모델 D)
           - 얇은 겔(50 μm) vs 두꺼운 겔(500 μm)
           - 두꺼운 겔에서 [E] 의존성 증가 → 확산 침투 제한
           - 얇은 겔에서 [E] 독립적 → 표면 반응 우세
        
        3. **생성물 첨가 실험** (모델 E)
           - 미리 절단된 펩타이드 조각 첨가
           - 반응 초기 속도 감소 → 생성물 억제 증명
           - 높은 [생성물]에서 α_max 감소 관찰
        
        4. **표면 처리 변화** (모델 F)
           - 양전하 표면 vs 음전하(PDA) vs 중성(PEG)
           - 음전하 표면에서 [E] 의존성 강화 → 흡착 증명
           - PEG 표면에서 흡착 감소 → 모델 D/E로 전환
        
        ### 🧬 고전적 메커니즘 테스트
        
        5. **Pulse-chase 실험** (모델 B)
           - t=5분에 신선한 효소 재투입
           - 곡선 재상승 → 기질 남음 (모델 A)
           - 변화 없음 → 효소 비활성화 (모델 B)
        
        6. **교반/유속 변화** (모델 C)
           - 정적 vs 회전 (100 rpm) vs 관류 (1 mL/min)
           - 유속 증가로 α 증가 → 물질전달 제한
           - 변화 없음 → 반응속도 제한 (모델 A/B)
        
        7. **기질 밀도 변화** (모델 A)
           - 0.5배, 1배, 2배 펩타이드 고정화
           - Fmax 비례 증가 → 기질 고갈
           - Fmax 불변 → 다른 메커니즘 우세
        
        8. **용액상 대조실험**
           - 가용성 기질 (같은 농도)
           - 완전 절단(α→1) → 표면/확산 문제
           - 불완전 절단 → 본질적 억제/비활성화
        
        ### 🎯 모델 결정 트리
        
        ```
        Fmax가 [E]에 따라 증가하는가?
        ├─ YES → 확장 모델 (D-F) 테스트
        │   ├─ 겔 두께 민감? → 모델 D (침투)
        │   ├─ 생성물 첨가로 감소? → 모델 E (억제)
        │   └─ 표면 전하 민감? → 모델 F (흡착)
        │
        └─ NO → 기본 모델 (A-C) 테스트
            ├─ Pulse-chase 반응? → 모델 A (기질)
            ├─ 시간에 따라 α_max↓? → 모델 B (비활성)
            └─ 유속에 민감? → 모델 C (확산)
        ```
        """)


