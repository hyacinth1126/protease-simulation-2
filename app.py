#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogel FRET Advanced Kinetic Analysis - Streamlit Application
"""

import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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

# Prep Raw Data 모드용 import
from prep_raw_data_mode.prep import (
    read_raw_data,
    fit_time_course,
    fit_calibration_curve,
    michaelis_menten_calibration,
    plot_calibration_curve
)
from data_interpolation_mode.interpolate_prism import (
    exponential_association,
    create_prism_interpolation_range
)
import os
import numpy as np
from scipy.optimize import curve_fit
import io
import base64

# Configure plotting
plt.rcParams['font.family'] = 'DejaVu Sans'
sns.set_style("whitegrid")


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


def main():
    """Main Streamlit app"""
    st.set_page_config(
        page_title="하이드로겔 FRET 고급 분석",
        page_icon="🔬",
        layout="wide"
    )
    
    st.title("🔬  Hydrogel FRET Simulation")
    st.markdown("---")
    
    # 모드 선택
    analysis_mode = st.sidebar.radio(
        "분석 모드 선택",
        ["Prep Raw Data 모드", "Data Interpolation 모드", "일반 분석 모드"],
        help="Prep Raw Data 모드: Raw data에서 MM Fitting | Data Interpolation 모드: Fitting 결과에서 Prism 스타일 보간 | 일반 분석 모드: 표준 FRET 분석"
    )
    
    st.markdown("---")
    
    # Prep Raw Data 모드
    if analysis_mode == "Prep Raw Data 모드":
        prep_raw_data_mode(st)
        return
    
    # Data Interpolation 모드
    if analysis_mode == "Data Interpolation 모드":
        data_interpolation_mode(st)
        return
    
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


if __name__ == "__main__":
    main()

