#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogel FRET Advanced Kinetic Analysis - Visualization Tools
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import List
from analysis import ModelResults


class Visualizer:
    """Visualization tools for model comparison"""
    
    @staticmethod
    def plot_raw_data(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_label: str = '시간 (초)'):
        """Plot raw fluorescence data with exponential fits and asymptotes"""
        fig = go.Figure()
        
        colors = px.colors.qualitative.Set1
        
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        for idx, conc in enumerate(sorted(df[conc_col].unique())):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            # Plot experimental data
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['FL_intensity'],
                mode='markers',
                name=f'{conc} {conc_unit} (데이터)',
                marker=dict(size=8, color=color),
                error_y=dict(type='data', array=subset['SD'], visible=True) if 'SD' in subset.columns else None
            ))
            
            # Plot exponential fit (if available)
            if 't_fit' in subset.columns and 'F_fit' in subset.columns:
                t_fit = subset['t_fit'].iloc[0]
                F_fit = subset['F_fit'].iloc[0]
                
                fig.add_trace(go.Scatter(
                    x=t_fit,
                    y=F_fit,
                    mode='lines',
                    name=f'{conc} {conc_unit} (지수 피팅)',
                    line=dict(color=color, dash='dash', width=2),
                    showlegend=True
                ))
                
                # Add asymptote line (Fmax)
                Fmax = subset['Fmax'].iloc[0]
                fig.add_trace(go.Scatter(
                    x=[t_fit[0], t_fit[-1]],
                    y=[Fmax, Fmax],
                    mode='lines',
                    name=f'{conc} {conc_unit} (점근선 Fmax)',
                    line=dict(color=color, dash='dot', width=1.5),
                    showlegend=True,
                    hovertemplate=f'Fmax = {Fmax:.1f}<extra></extra>'
                ))
        
        fig.update_layout(
            title='원본 형광 데이터 with 지수 피팅 및 점근선',
            xaxis_title=time_label,
            yaxis_title='형광 강도 (RFU)',
            template='plotly_white',
            hovermode='x unified',
            height=600
        )
        
        return fig
    
    @staticmethod
    def plot_normalized_data(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_label: str = '시간 (초)'):
        """Plot normalized data (fraction cleaved)"""
        fig = go.Figure()
        
        colors = px.colors.qualitative.Set1
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        concentrations = sorted(df[conc_col].unique())
        
        for idx, conc in enumerate(concentrations):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            # Plot data points
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['alpha'],
                mode='markers',
                name=f'{conc} {conc_unit}',
                marker=dict(size=8, color=color),
                legendgroup=f'group{idx}'
            ))
        
        fig.update_layout(
            title='정규화 데이터: 절단 비율 α(t)',
            xaxis_title=time_label,
            yaxis_title='절단 비율 α',
            template='plotly_white',
            hovermode='x unified',
            height=600
        )
        
        return fig
    
    @staticmethod
    def plot_model_fits(df: pd.DataFrame, results: List[ModelResults], 
                       conc_unit: str = 'μg/mL', time_label: str = '시간 (초)'):
        """Plot all model fits together"""
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('모델 피팅', '잔차'),
            vertical_spacing=0.15,
            row_heights=[0.7, 0.3]
        )
        
        colors = px.colors.qualitative.Set1
        
        # Plot experimental data
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        for idx, conc in enumerate(sorted(df[conc_col].unique())):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['alpha'],
                mode='markers',
                name=f'{conc} {conc_unit} (데이터)',
                marker=dict(size=6, color=color),
                showlegend=(idx == 0)
            ), row=1, col=1)
        
        # Plot model predictions
        line_styles = ['solid', 'dash', 'dot']
        for model_idx, result in enumerate(results):
            if result is None:
                continue
            
            for idx, conc in enumerate(sorted(df['enzyme_ugml'].unique())):
                subset = df[df['enzyme_ugml'] == conc]
                indices = subset.index
                
                fig.add_trace(go.Scatter(
                    x=subset['time_s'],
                    y=result.predictions[indices],
                    mode='lines',
                    name=result.name if idx == 0 else None,
                    line=dict(width=2, dash=line_styles[model_idx % 3]),
                    showlegend=(idx == 0),
                    legendgroup=result.name
                ), row=1, col=1)
        
        # Plot residuals
        for model_idx, result in enumerate(results):
            if result is None:
                continue
            
            fig.add_trace(go.Scatter(
                x=df['time_s'],
                y=result.residuals,
                mode='markers',
                name=f'{result.name} 잔차',
                marker=dict(size=4),
                showlegend=False
            ), row=2, col=1)
        
        fig.add_hline(y=0, line_dash="dash", line_color="gray", row=2, col=1)
        
        fig.update_xaxes(title_text=time_label, row=2, col=1)
        fig.update_yaxes(title_text="α", row=1, col=1)
        fig.update_yaxes(title_text="잔차", row=2, col=1)
        
        fig.update_layout(
            height=800,
            template='plotly_white',
            hovermode='x unified'
        )
        
        return fig
    
    @staticmethod
    def plot_initial_rates(df: pd.DataFrame, conc_unit: str = 'μg/mL', time_unit: str = 's'):
        """Plot initial rates v0 vs [E] to check linearity"""
        # Calculate initial slopes (0-2 time units)
        cutoff_time = 2 if time_unit == 's' else 0.5  # 2 seconds or 0.5 minutes
        initial_data = df[df['time_s'] <= cutoff_time].copy()
        
        rates = []
        concentrations = []
        
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        for conc in sorted(df[conc_col].unique()):
            subset = initial_data[initial_data[conc_col] == conc]
            if len(subset) >= 2:
                # Linear fit to get slope
                coeffs = np.polyfit(subset['time_s'], subset['alpha'], 1)
                v0 = coeffs[0]  # slope = dα/dt
                rates.append(v0)
                concentrations.append(conc)
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=concentrations,
            y=rates,
            mode='markers+lines',
            name='초기 속도',
            marker=dict(size=10)
        ))
        
        # Linear fit
        if len(concentrations) >= 2:
            coeffs = np.polyfit(concentrations, rates, 1)
            x_fit = np.linspace(min(concentrations), max(concentrations), 100)
            y_fit = coeffs[0] * x_fit + coeffs[1]
            
            fig.add_trace(go.Scatter(
                x=x_fit,
                y=y_fit,
                mode='lines',
                name=f'선형 피팅 (R² = {np.corrcoef(concentrations, rates)[0,1]**2:.4f})',
                line=dict(dash='dash')
            ))
        
        time_unit_label = "초" if time_unit == 's' else "분"
        time_unit_abbr = "s⁻¹" if time_unit == 's' else "min⁻¹"
        fig.update_layout(
            title=f'초기 속도 분석 (0-{cutoff_time}{time_unit_label})',
            xaxis_title=f'[효소] ({conc_unit})',
            yaxis_title=f'v₀ ({time_unit_abbr})',
            template='plotly_white'
        )
        
        return fig
    
    @staticmethod
    def plot_individual_model(df: pd.DataFrame, result: ModelResults, 
                             conc_unit: str = 'μg/mL', time_label: str = '시간 (초)',
                             model_color: str = '#FF6B6B'):
        """Plot a single model fit with experimental data"""
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=(f'{result.name} 피팅 결과', '잔차'),
            vertical_spacing=0.15,
            row_heights=[0.7, 0.3]
        )
        
        colors = px.colors.qualitative.Set1
        
        # Plot experimental data
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        for idx, conc in enumerate(sorted(df[conc_col].unique())):
            subset = df[df[conc_col] == conc]
            color = colors[idx % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=subset['alpha'],
                mode='markers',
                name=f'{conc} {conc_unit}',
                marker=dict(size=8, color=color),
                legendgroup=f'conc_{idx}'
            ), row=1, col=1)
            
            # Plot model fit for this concentration
            indices = subset.index
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=result.predictions[indices],
                mode='lines',
                name=f'{conc} {conc_unit} (피팅)',
                line=dict(width=3, color=color),
                showlegend=False,
                legendgroup=f'conc_{idx}'
            ), row=1, col=1)
            
            # Plot residuals for this concentration
            fig.add_trace(go.Scatter(
                x=subset['time_s'],
                y=result.residuals[indices],
                mode='markers',
                marker=dict(size=5, color=color),
                showlegend=False,
                legendgroup=f'conc_{idx}'
            ), row=2, col=1)
        
        fig.add_hline(y=0, line_dash="dash", line_color="gray", row=2, col=1)
        
        fig.update_xaxes(title_text=time_label, row=2, col=1)
        fig.update_yaxes(title_text="α", row=1, col=1)
        fig.update_yaxes(title_text="잔차", row=2, col=1)
        
        # Add statistics annotation
        stats_text = f"<b>피팅 통계</b><br>R² = {result.r_squared:.4f}<br>RMSE = {result.rmse:.4f}"
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.02, y=0.98,
            xanchor='left', yanchor='top',
            text=stats_text,
            showarrow=False,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor=model_color,
            borderwidth=2,
            font=dict(size=11),
            row=1, col=1
        )
        
        fig.update_layout(
            height=800,
            template='plotly_white',
            hovermode='x unified',
            showlegend=True
        )
        
        return fig
    
    @staticmethod
    def create_comparison_table(results: List[ModelResults]) -> pd.DataFrame:
        """Create comparison table for models"""
        data = []
        
        for result in results:
            if result is None:
                continue
            
            row = {
                '모델': result.name,
                'R²': f"{result.r_squared:.4f}",
                'RMSE': f"{result.rmse:.4f}",
                'AIC': f"{result.aic:.2f}",
                'BIC': f"{result.bic:.2f}",
                '파라미터 수': len(result.params)
            }
            
            # Add parameter values
            for param_name, param_value in result.params.items():
                if param_name == 'kcat_KM':
                    row[param_name] = f"{param_value:.2e} M⁻¹s⁻¹"
                elif param_name == 'kd':
                    row[param_name] = f"{param_value:.4f} s⁻¹"
                elif param_name == 'km':
                    row[param_name] = f"{param_value:.2e} m/s"
                elif param_name == 'Gamma0':
                    row[param_name] = f"{param_value:.2f} pmol/cm²"
                elif param_name == 'alpha_inf':
                    row[param_name] = f"{param_value:.4f}"
                elif param_name == 'k_access':
                    row[param_name] = f"{param_value:.2e} M⁻¹"
                elif param_name == 'Ki_eff':
                    row[param_name] = f"{param_value:.4f}"
                elif param_name == 'k_ads':
                    row[param_name] = f"{param_value:.4f} s⁻¹"
                elif param_name == 'K_ads':
                    row[param_name] = f"{param_value:.2e} M⁻¹"
                else:
                    row[param_name] = f"{param_value:.4e}"
            
            data.append(row)
        
        return pd.DataFrame(data)

