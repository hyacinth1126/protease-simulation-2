#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# .\.venv\Scripts\python.exe -m streamlit run app.py .\.venv\Scripts\python.exe -m streamlit run app.py
# author: hyacinth1126
"""
Hydrogel FRET Advanced Kinetic Analysis - Streamlit Application
"""

import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns

# UI 모듈 import
from app_ui.prep_mode import prep_raw_data_mode
from app_ui.interp_mode import data_interpolation_mode
from app_ui.general_analysis_mode import general_analysis_mode
from app_ui.footer import render_footer

# Configure plotting
plt.rcParams['font.family'] = 'DejaVu Sans'
sns.set_style("whitegrid")


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
    # 항상 하단에 푸터 렌더링
    render_footer()
    
    # Prep Raw Data 모드
    if analysis_mode == "Prep Raw Data 모드":
        prep_raw_data_mode(st)
        return
    
    # Data Interpolation 모드
    if analysis_mode == "Data Interpolation 모드":
        data_interpolation_mode(st)
        return
    
    # 일반 분석 모드
    general_analysis_mode(st)


if __name__ == "__main__":
    main()
