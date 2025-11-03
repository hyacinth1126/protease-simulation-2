#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogel FRET Advanced Kinetic Analysis - Core Logic
Three competing models: Substrate Depletion, Enzyme Deactivation, Mass-Transfer Limitation
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from dataclasses import dataclass
from typing import Dict
import warnings

# Suppress scipy optimization warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, module='scipy')
warnings.filterwarnings('ignore', message='Covariance of the parameters could not be estimated')


@dataclass
class ModelResults:
    """Container for model fitting results"""
    name: str
    params: Dict[str, float]
    params_std: Dict[str, float]
    aic: float
    bic: float
    r_squared: float
    rmse: float
    predictions: np.ndarray
    residuals: np.ndarray


class UnitStandardizer:
    """
    Step 2: Standardize units
    
    1) Time: min/sec → time_s (seconds)
    2) Concentration: μg/mL, ng/mL → μM, nM (molar concentration)
       - Requires molecular weight (MW) for conversion
    3) Fluorescence: RFU or FL_intensity (no conversion, just column selection)
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        """
        Parameters:
        - enzyme_mw: Molecular weight in kDa (default: 56.6 kDa for Kgp)
        """
        self.enzyme_mw = enzyme_mw  # kDa
    
    def standardize(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Standardize units in the dataframe
        
        Steps:
        1. Time standardization: min → s
        2. Concentration standardization: mass → molar (if needed)
        3. Fluorescence column selection: RFU or FL_intensity
        
        Parameters:
        - df: Raw DataFrame
        
        Returns:
        - DataFrame with standardized units
        """
        df_std = df.copy()
        
        # Step 1: Time standardization
        if 'time_min' in df_std.columns:
            df_std['time_s'] = df_std['time_min'] * 60  # min → s
        elif 'time_s' not in df_std.columns:
            # Try to find time column
            time_cols = [col for col in df_std.columns if 'time' in col.lower()]
            if len(time_cols) > 0:
                df_std['time_s'] = df_std[time_cols[0]]
        
        # Step 2: Fluorescence column standardization
        if 'RFU' in df_std.columns and 'FL_intensity' not in df_std.columns:
            df_std['FL_intensity'] = df_std['RFU']
        elif 'FL_intensity' not in df_std.columns:
            # Try to find fluorescence column
            fluor_cols = [col for col in df_std.columns 
                         if any(keyword in col.upper() for keyword in ['FLUOR', 'RFU', 'FL_', 'INTENSITY'])]
            if len(fluor_cols) > 0:
                df_std['FL_intensity'] = df_std[fluor_cols[0]]
        
        # Step 3: Concentration standardization (mass → molar)
        conc_col = None
        conc_unit_type = None
        
        # Check for existing molar concentration columns
        if 'peptide_uM' in df_std.columns:
            conc_col = 'peptide_uM'
            conc_unit_type = 'molar'
            target_col = 'enzyme_uM'
        elif 'enzyme_uM' in df_std.columns:
            conc_col = 'enzyme_uM'
            conc_unit_type = 'molar'
            target_col = 'enzyme_uM'
        elif 'E_nM' in df_std.columns:
            conc_col = 'E_nM'
            conc_unit_type = 'molar'
            target_col = 'E_nM'
        # Check for mass concentration columns
        elif 'enzyme_ugml' in df_std.columns or 'enzyme_ug/mL' in df_std.columns:
            conc_col = 'enzyme_ugml' if 'enzyme_ugml' in df_std.columns else 'enzyme_ug/mL'
            conc_unit_type = 'mass'
            target_col = 'enzyme_nM'  # Changed from uM to nM for μg/mL range
        elif 'peptide_ugml' in df_std.columns or 'peptide_ug/mL' in df_std.columns:
            conc_col = 'peptide_ugml' if 'peptide_ugml' in df_std.columns else 'peptide_ug/mL'
            conc_unit_type = 'mass'
            target_col = 'peptide_nM'  # Changed from uM to nM for μg/mL range
        elif 'enzyme_ngml' in df_std.columns or 'enzyme_ng/mL' in df_std.columns:
            conc_col = 'enzyme_ngml' if 'enzyme_ngml' in df_std.columns else 'enzyme_ng/mL'
            conc_unit_type = 'mass'
            target_col = 'enzyme_pM'  # Changed from nM to pM for ng/mL range
        
        # Convert mass to molar if needed
        if conc_col and conc_unit_type == 'mass':
            MW_kDa = self.enzyme_mw  # kDa
            
            if 'uM' in target_col:
                # Convert μg/mL to μM
                # Formula: μM = (μg/mL) / (MW in kDa) / 1000
                df_std[target_col] = df_std[conc_col] / MW_kDa / 1000
            elif 'nM' in target_col:
                # Convert μg/mL to nM
                # Formula: nM = (μg/mL) / (MW in kDa)
                # Example: 0.3125 μg/mL / 56.6 kDa = 5.52 nM
                df_std[target_col] = df_std[conc_col] / MW_kDa
            elif 'pM' in target_col:
                # Convert ng/mL to pM
                # Formula: pM = (ng/mL) / (MW in kDa)
                # Example: 3.125 ng/mL / 56.6 kDa = 55.2 pM
                df_std[target_col] = df_std[conc_col] / MW_kDa
            
            # Store original column name for reference
            df_std['_original_conc_col'] = conc_col
        elif conc_col and conc_unit_type == 'molar':
            # Already in molar, just ensure target_col exists
            if target_col not in df_std.columns:
                df_std[target_col] = df_std[conc_col]
        
        return df_std


class RegionDivider:
    """
    Step 4: Divide kinetic regions
    
    Divides normalized data into:
    1. Initial linear region
    2. Exponential growth region  
    3. Plateau region
    
    Note: Specific logic to be provided later by user
    """
    
    def __init__(self):
        pass
    
    def divide_regions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Divide kinetic data into regions
        
        Parameters:
        - df: Normalized DataFrame (must have 'alpha' or 'alpha_temp', 'time_s', etc.)
        
        Returns:
        - DataFrame with additional 'region' column
        """
        df_regions = df.copy()
        
        # Use alpha if available, otherwise use alpha_temp
        alpha_col = 'alpha' if 'alpha' in df_regions.columns else 'alpha_temp'
        
        # TODO: Implement region division logic
        # For now, just add placeholder column
        df_regions['region'] = 'unknown'
        
        return df_regions


class DataNormalizer:
    """
    Step 3: Normalization (two stages)
    
    Stage 1: Temporary normalization (model-free threshold)
    - F0 = minimum fluorescence value
    - Fmax = maximum fluorescence value
    
    Stage 2: Final normalization (after region division)
    - F0 = minimum fluorescence value (same as stage 1)
    - Fmax = plateau region average (or F∞ from exponential region if no plateau)
    """
    
    def __init__(self):
        self.time_col = None
        self.conc_col = None
        self.fluor_col = None
    
    def _detect_columns(self, df: pd.DataFrame):
        """Detect and store column names"""
        # Time column
        self.time_col = 'time_s' if 'time_s' in df.columns else 'time_min'
        
        # Concentration column - prioritize molar concentration units
        if 'peptide_uM' in df.columns:
            self.conc_col = 'peptide_uM'
            self.conc_unit_type = 'molar'  # uM
        elif 'enzyme_uM' in df.columns:
            self.conc_col = 'enzyme_uM'
            self.conc_unit_type = 'molar'  # uM
        elif 'E_nM' in df.columns:
            self.conc_col = 'E_nM'
            self.conc_unit_type = 'molar'  # nM
        elif 'enzyme_ugml' in df.columns:
            self.conc_col = 'enzyme_ugml'
            self.conc_unit_type = 'mass'  # μg/mL
        else:
            # Try to find any column that might be concentration
            conc_cols = [col for col in df.columns 
                        if 'time' not in col.lower() and 'fluor' not in col.lower() 
                        and 'FL_' not in col and 'intensity' not in col.lower() 
                        and 'RFU' not in col and 'region' not in col.lower()]
            if len(conc_cols) > 0:
                self.conc_col = conc_cols[0]
                self.conc_unit_type = 'unknown'
            else:
                self.conc_col = 'concentration'
                self.conc_unit_type = 'unknown'
        
        # Fluorescence column
        self.fluor_col = 'FL_intensity' if 'FL_intensity' in df.columns else 'RFU' if 'RFU' in df.columns else 'fluor'
    
    def normalize_temporary(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Stage 1: Temporary normalization using model-free threshold method
        
        F0 = minimum fluorescence value (per concentration)
        Fmax = maximum fluorescence value (per concentration)
        
        Parameters:
        - df: DataFrame with standardized units
        
        Returns:
        - DataFrame with temporary alpha, F0_temp, Fmax_temp columns
        """
        self._detect_columns(df)
        df_temp = df.copy()
        
        def temp_normalize_group(g):
            g_sorted = g.sort_values(self.time_col)
            
            # F0: minimum fluorescence value
            F0_temp = float(g_sorted[self.fluor_col].min())
            
            # Fmax: maximum fluorescence value
            Fmax_temp = float(g_sorted[self.fluor_col].max())
            
            # Ensure Fmax > F0
            if Fmax_temp <= F0_temp:
                Fmax_temp = F0_temp + 100
            
            # Calculate temporary alpha
            g = g_sorted.copy()
            g['alpha_temp'] = (g[self.fluor_col] - F0_temp) / (Fmax_temp - F0_temp)
            g['alpha_temp'] = g['alpha_temp'].clip(0, 1)
            g['F0_temp'] = F0_temp
            g['Fmax_temp'] = Fmax_temp
            
            # Store column names
            g['time_col_name'] = self.time_col
            g['conc_col_name'] = self.conc_col
            g['fluor_col_name'] = self.fluor_col
            
            return g
        
        df_normalized = df_temp.groupby(self.conc_col, group_keys=False).apply(temp_normalize_group)
        return df_normalized
    
    def normalize_final(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Stage 2: Final normalization using region information
        
        F0 = minimum fluorescence value (same as temporary)
        Fmax = plateau region average (or F∞ from exponential region if no plateau)
        
        Parameters:
        - df: DataFrame with temporary normalization and region division
        
        Returns:
        - DataFrame with final alpha, F0, Fmax columns
        """
        self._detect_columns(df)
        df_final = df.copy()
        
        def final_normalize_group(g):
            g_sorted = g.sort_values(self.time_col)
            
            # F0: same as temporary (minimum)
            F0 = float(g_sorted['F0_temp'].iloc[0])
            
            # Fmax: determine based on region
            if 'region' not in g_sorted.columns:
                # No region info, fallback to temporary Fmax
                Fmax = float(g_sorted['Fmax_temp'].iloc[0])
                region_used = 'fallback_temp'
            else:
                # Check if plateau region exists
                plateau_data = g_sorted[g_sorted['region'] == 'plateau']
                
                if len(plateau_data) > 0:
                    # Use plateau region average
                    Fmax = float(plateau_data[self.fluor_col].mean())
                    region_used = 'plateau_mean'
                else:
                    # No plateau, fit exponential to exponential region
                    exp_data = g_sorted[g_sorted['region'] == 'exponential']
                    
                    if len(exp_data) >= 3:
                        # Fit exponential to determine F∞
                        try:
                            t_data = exp_data[self.time_col].values
                            F_data = exp_data[self.fluor_col].values
                            
                            def exp_func(t, A, k):
                                return F0 + A * (1 - np.exp(-k * t))
                            
                            A_guess = F_data[-1] - F0 if len(F_data) > 0 else 1000
                            p0 = [A_guess, 0.1]
                            
                            popt, _ = curve_fit(
                                exp_func, t_data, F_data,
                                p0=p0,
                                bounds=([0, 0.001], [np.inf, 10]),
                                maxfev=5000
                            )
                            
                            A_fit, _ = popt
                            Fmax = F0 + A_fit  # F∞
                            region_used = 'exponential_Finf'
                        except:
                            # Fallback to maximum
                            Fmax = float(g_sorted[self.fluor_col].max())
                            region_used = 'fallback_max'
                    else:
                        # Not enough exponential data, use maximum
                        Fmax = float(g_sorted[self.fluor_col].max())
                        region_used = 'fallback_max'
            
            # Ensure Fmax > F0
            if Fmax <= F0:
                Fmax = F0 + 100
            
            # Calculate final alpha
            g = g_sorted.copy()
            g['alpha'] = (g[self.fluor_col] - F0) / (Fmax - F0)
            g['alpha'] = g['alpha'].clip(0, 1)
            g['F0'] = F0
            g['Fmax'] = Fmax
            g['Fmax_method'] = region_used
            
            return g
        
        df_normalized = df_final.groupby(self.conc_col, group_keys=False).apply(final_normalize_group)
        return df_normalized


class ModelA_SubstrateDepletion:
    """
    Model A: Substrate Depletion
    
    Surface reaction with immobilized substrate:
    dΓ/dt = -v = -(kcat * Es * Γ) / (KM + Γ)
    
    Simplified (Γ << KM): dΓ/dt = -kobs * Γ
    where kobs ≈ (kcat/KM) * Es
    
    Solution: α(t) = 1 - exp(-kobs * t)
    With mass transfer: Es ≈ Eb / (1 + Da), Da = (kcat*Γ0)/(KM*km)
    
    Parameters to fit:
    - kcat_KM: catalytic efficiency (M^-1 s^-1)
    - Gamma0: initial surface substrate density (pmol/cm²)
    - km: mass transfer coefficient (optional, m/s)
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        self.enzyme_mw = enzyme_mw  # kDa
        self.name = "Model A: Substrate Depletion"
        
    def model_simple(self, t: np.ndarray, E_M: float, kcat_KM: float) -> np.ndarray:
        """
        Simplified first-order kinetics
        α(t) = 1 - exp(-kobs * t)
        kobs = (kcat/KM) * [E]
        """
        kobs = kcat_KM * E_M
        alpha = 1.0 - np.exp(-kobs * t)
        return alpha
    
    def model_with_saturation(self, t: np.ndarray, E_M: float, 
                             kcat_KM: float, Gamma0: float, KM: float = 1e-6) -> np.ndarray:
        """
        Full Michaelis-Menten kinetics
        Numerically integrate: dΓ/dt = -(kcat/KM) * E * Γ
        """
        kcat = kcat_KM * KM
        
        def dydt(t_val, y):
            Gamma = y[0]
            if Gamma < 0:
                return [0]
            v = (kcat * E_M * Gamma) / (KM + Gamma)
            return [-v]
        
        sol = solve_ivp(dydt, [t[0], t[-1]], [Gamma0], t_eval=t, method='LSODA')
        Gamma_t = sol.y[0]
        alpha = 1.0 - Gamma_t / Gamma0
        return np.clip(alpha, 0, 1)
    
    def fit_global(self, df: pd.DataFrame, verbose_callback=None) -> ModelResults:
        """
        Global fit to all concentration data
        
        Parameters:
        - df: DataFrame with normalized data
        - verbose_callback: Optional callback function(message) for logging
        
        Returns:
        - ModelResults object or None if fitting fails
        """
        # Use all data - no filtering
        df_fit = df.copy()
        
        if verbose_callback:
            verbose_callback(f"Model A: 전체 {len(df_fit)}개 데이터 포인트 사용 (필터링 없음)")
        
        if len(df_fit) < 5:
            if verbose_callback:
                verbose_callback(f"Model A: 피팅 데이터가 충분하지 않습니다. (사용 가능: {len(df_fit)}개)", level="error")
            return None
        
        # Convert enzyme concentration to M
        MW = self.enzyme_mw * 1000  # g/mol
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        
        if 'uM' in conc_col:
            # Already in μM, convert to M
            df_fit['E_M'] = df_fit[conc_col] * 1e-6
        elif 'nM' in conc_col:
            # Already in nM, convert to M
            df_fit['E_M'] = df_fit[conc_col] * 1e-9
        else:
            # Mass concentration (μg/mL), convert to M
            df_fit['E_M'] = (df_fit[conc_col] / MW) * 1e-6
        
        t = df_fit['time_s'].values
        E = df_fit['E_M'].values
        y = df_fit['alpha'].values
        
        # Store globally for model function access
        self._t_data = t
        self._E_data = E
        
        # Fit simplified model: α = 1 - exp(-(kcat/KM)*E*t)
        def model_func(x_dummy, kcat_KM):
            """Model function for curve_fit"""
            result = np.zeros_like(self._t_data)
            for i in range(len(self._t_data)):
                kobs = kcat_KM * self._E_data[i]
                result[i] = 1.0 - np.exp(-kobs * self._t_data[i])
            return result
        
        try:
            # Better initial guess estimation
            kobs_estimates = []
            for conc in df_fit['enzyme_ugml'].unique():
                subset = df_fit[df_fit['enzyme_ugml'] == conc].sort_values('time_s')
                if len(subset) >= 3:
                    early = subset.head(min(5, len(subset)))
                    t_early = early['time_s'].values
                    alpha_early = early['alpha'].values
                    
                    alpha_safe = np.clip(alpha_early, 0.01, 0.99)
                    y_log = -np.log(1 - alpha_safe)
                    
                    if len(t_early) >= 2 and t_early[-1] > t_early[0]:
                        kobs_est = (y_log[-1] - y_log[0]) / (t_early[-1] - t_early[0])
                        E_conc = subset['E_M'].iloc[0]
                        if E_conc > 0 and kobs_est > 0:
                            kcat_KM_est = kobs_est / E_conc
                            kobs_estimates.append(kcat_KM_est)
            
            if len(kobs_estimates) > 0:
                kcat_KM_guess = np.median(kobs_estimates)
                kcat_KM_guess = np.clip(kcat_KM_guess, 1e2, 1e10)
            else:
                kcat_KM_guess = 1e5
            
            p0 = [kcat_KM_guess]
            
            if verbose_callback:
                verbose_callback(f"Model A 초기값: kcat/KM = {kcat_KM_guess:.2e} M⁻¹s⁻¹")
            
            # Dummy x values (just indices)
            x_dummy = np.arange(len(y))
            
            # Curve fit with wider bounds
            popt, pcov = curve_fit(
                model_func, x_dummy, y, 
                p0=p0, bounds=([1e3], [1e9]), maxfev=20000
            )
            
            kcat_KM_fit = popt[0]
            perr = np.sqrt(np.diag(pcov))
            kcat_KM_std = perr[0]
            
            # Calculate predictions and metrics
            y_pred = model_func(x_dummy, kcat_KM_fit)
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(np.mean(residuals ** 2))
            
            # AIC/BIC
            n = len(y)
            k = 1  # number of parameters
            if ss_res > 0:
                log_likelihood = -n/2 * np.log(2*np.pi*ss_res/n) - n/2
            else:
                log_likelihood = 0
            aic = 2*k - 2*log_likelihood
            bic = k*np.log(n) - 2*log_likelihood
            
            params = {'kcat_KM': kcat_KM_fit}
            params_std = {'kcat_KM': kcat_KM_std}
            
            # Generate full predictions for all data (including saturated)
            t_full = df['time_s'].values
            conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
            
            if 'uM' in conc_col:
                E_full = df[conc_col].values * 1e-6
            elif 'nM' in conc_col:
                E_full = df[conc_col].values * 1e-9
            else:
                E_full = (df[conc_col].values / MW) * 1e-6
            y_pred_full = np.zeros(len(t_full))
            for i in range(len(t_full)):
                kobs = kcat_KM_fit * E_full[i]
                y_pred_full[i] = 1.0 - np.exp(-kobs * t_full[i])
            
            return ModelResults(
                name=self.name,
                params=params,
                params_std=params_std,
                aic=aic,
                bic=bic,
                r_squared=r_squared,
                rmse=rmse,
                predictions=y_pred_full,
                residuals=df['alpha'].values - y_pred_full
            )
            
        except Exception as e:
            if verbose_callback:
                verbose_callback(f"Model A 피팅 오류: {str(e)}", level="error")
                import traceback
                verbose_callback(traceback.format_exc(), level="debug")
            return None


class ModelB_EnzymeDeactivation:
    """
    Model B: Enzyme Deactivation
    
    Active enzyme decays: dE/dt = -kd * E → E(t) = E0 * exp(-kd*t)
    Surface reaction: dΓ/dt = -(kcat/KM) * E(t) * Γ
    
    For Γ << KM (first order in Γ):
    α(t) = α∞ * [1 - exp(-kobs*t)], where α∞ < 1 if enzyme deactivates fast
    
    Or integrated: α(t) = (kcat/KM)*E0/kd * [1 - exp(-kd*t)] when Γ>>KM
    
    Approximate: α(t) ≈ 1 - exp[-(kcat/KM)*E0/kd * (1 - exp(-kd*t))]
    
    Parameters to fit:
    - kcat_KM: catalytic efficiency (M^-1 s^-1)
    - kd: enzyme deactivation rate (s^-1)
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        self.enzyme_mw = enzyme_mw
        self.name = "Model B: Enzyme Deactivation"
    
    def model(self, t: np.ndarray, E0_M: float, kcat_KM: float, kd: float) -> np.ndarray:
        """
        Model with enzyme deactivation
        Assuming first-order in both enzyme and substrate
        
        Solution: α(t) = 1 - exp[-(kcat/KM)*E0/kd * (1 - exp(-kd*t))]
        """
        if kd <= 0:
            # No deactivation, revert to Model A
            kobs = kcat_KM * E0_M
            return 1.0 - np.exp(-kobs * t)
        
        # Integrated form
        integral_term = (kcat_KM * E0_M / kd) * (1 - np.exp(-kd * t))
        alpha = 1.0 - np.exp(-integral_term)
        return np.clip(alpha, 0, 1)
    
    def fit_global(self, df: pd.DataFrame, verbose_callback=None) -> ModelResults:
        """Global fit including enzyme deactivation"""
        df_fit = df.copy()
        
        if verbose_callback:
            verbose_callback(f"Model B: 전체 {len(df_fit)}개 데이터 포인트 사용 (필터링 없음)")
        
        if len(df_fit) < 5:
            if verbose_callback:
                verbose_callback(f"Model B: 피팅 데이터가 충분하지 않습니다. (사용 가능: {len(df_fit)}개)", level="error")
            return None
        
        MW = self.enzyme_mw * 1000
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        
        if 'uM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-6
        elif 'nM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-9
        else:
            df_fit['E_M'] = (df_fit[conc_col] / MW) * 1e-6
        
        t = df_fit['time_s'].values
        E = df_fit['E_M'].values
        y = df_fit['alpha'].values
        
        self._t_data = t
        self._E_data = E
        
        def model_func(x_dummy, kcat_KM, kd):
            """Model function for curve_fit"""
            result = np.zeros_like(self._t_data)
            for i in range(len(self._t_data)):
                result[i] = self.model(self._t_data[i], self._E_data[i], kcat_KM, kd)
            return result
        
        try:
            # Better initial guess estimation
            kobs_estimates = []
            for conc in df_fit['enzyme_ugml'].unique():
                subset = df_fit[df_fit['enzyme_ugml'] == conc].sort_values('time_s')
                if len(subset) >= 3:
                    early = subset.head(min(5, len(subset)))
                    t_early = early['time_s'].values
                    alpha_early = early['alpha'].values
                    alpha_safe = np.clip(alpha_early, 0.01, 0.99)
                    y_log = -np.log(1 - alpha_safe)
                    if len(t_early) >= 2 and t_early[-1] > t_early[0]:
                        kobs_est = (y_log[-1] - y_log[0]) / (t_early[-1] - t_early[0])
                        E_conc = subset['E_M'].iloc[0]
                        if E_conc > 0 and kobs_est > 0:
                            kcat_KM_est = kobs_est / E_conc
                            kobs_estimates.append(kcat_KM_est)
            
            if len(kobs_estimates) > 0:
                kcat_KM_guess = np.median(kobs_estimates)
                kcat_KM_guess = np.clip(kcat_KM_guess, 1e2, 1e10)
            else:
                kcat_KM_guess = 1e5
            
            p0 = [kcat_KM_guess, 0.5]
            
            if verbose_callback:
                verbose_callback(f"Model B 초기값: kcat/KM = {kcat_KM_guess:.2e} M⁻¹s⁻¹, kd = 0.5 s⁻¹")
            
            x_dummy = np.arange(len(y))
            
            popt, pcov = curve_fit(
                model_func, x_dummy, y, p0=p0, 
                bounds=([1e3, 0.001], [1e9, 100]),
                maxfev=20000
            )
            
            kcat_KM_fit, kd_fit = popt
            perr = np.sqrt(np.diag(pcov))
            
            y_pred = model_func(x_dummy, kcat_KM_fit, kd_fit)
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(np.mean(residuals ** 2))
            
            n = len(y)
            k = 2
            if ss_res > 0:
                log_likelihood = -n/2 * np.log(2*np.pi*ss_res/n) - n/2
            else:
                log_likelihood = 0
            aic = 2*k - 2*log_likelihood
            bic = k*np.log(n) - 2*log_likelihood
            
            params = {'kcat_KM': kcat_KM_fit, 'kd': kd_fit}
            params_std = {'kcat_KM': perr[0], 'kd': perr[1]}
            
            t_full = df['time_s'].values
            conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
            
            if 'uM' in conc_col:
                E_full = df[conc_col].values * 1e-6
            elif 'nM' in conc_col:
                E_full = df[conc_col].values * 1e-9
            else:
                E_full = (df[conc_col].values / MW) * 1e-6
            y_pred_full = np.zeros(len(t_full))
            for i in range(len(t_full)):
                y_pred_full[i] = self.model(t_full[i], E_full[i], kcat_KM_fit, kd_fit)
            
            return ModelResults(
                name=self.name,
                params=params,
                params_std=params_std,
                aic=aic,
                bic=bic,
                r_squared=r_squared,
                rmse=rmse,
                predictions=y_pred_full,
                residuals=df['alpha'].values - y_pred_full
            )
            
        except Exception as e:
            if verbose_callback:
                verbose_callback(f"Model B 피팅 오류: {str(e)}", level="error")
                import traceback
                verbose_callback(traceback.format_exc(), level="debug")
            return None


class ModelC_MassTransfer:
    """
    Model C: Mass-Transfer Limitation
    
    Diffusion-limited enzyme delivery to surface
    Simplified using mass transfer coefficient km:
    
    Flux: J = km * (Eb - Es)
    Surface reaction: v = (kcat/KM) * Es * Γ
    
    Damköhler number: Da = (kcat*Γ0) / (KM*km)
    
    For high Da (reaction >> diffusion):
    Es << Eb, reaction is diffusion-controlled
    
    Approximate model: α(t) = 1 - exp(-keff*t)
    where keff = (kcat/KM)*Es, and Es depends on km and Da
    
    Parameters to fit:
    - kcat_KM: catalytic efficiency
    - km: mass transfer coefficient (m/s)
    - Gamma0: surface substrate density
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        self.enzyme_mw = enzyme_mw
        self.name = "Model C: Mass-Transfer Limitation"
    
    def calculate_Es(self, Eb: float, kcat_KM: float, Gamma0: float, 
                     km: float, KM: float = 1e-6) -> float:
        """
        Calculate surface enzyme concentration accounting for mass transfer
        Es ≈ Eb / (1 + Da), where Da = (kcat*Γ0)/(KM*km)
        """
        kcat = kcat_KM * KM
        km_cm = km * 100
        Da = (kcat * Gamma0 * 1e-12) / (KM * km_cm)
        Es = Eb / (1 + Da)
        return Es
    
    def model(self, t: np.ndarray, Eb_M: float, kcat_KM: float, 
              km: float, Gamma0: float = 1.0) -> np.ndarray:
        """Model with mass transfer limitation"""
        Es = self.calculate_Es(Eb_M, kcat_KM, Gamma0, km)
        kobs = kcat_KM * Es
        alpha = 1.0 - np.exp(-kobs * t)
        return np.clip(alpha, 0, 1)
    
    def fit_global(self, df: pd.DataFrame, verbose_callback=None) -> ModelResults:
        """Global fit with mass transfer"""
        df_fit = df.copy()
        
        if verbose_callback:
            verbose_callback(f"Model C: 전체 {len(df_fit)}개 데이터 포인트 사용 (필터링 없음)")
        
        if len(df_fit) < 5:
            if verbose_callback:
                verbose_callback(f"Model C: 피팅 데이터가 충분하지 않습니다. (사용 가능: {len(df_fit)}개)", level="error")
            return None
        
        MW = self.enzyme_mw * 1000
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        
        if 'uM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-6
        elif 'nM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-9
        else:
            df_fit['E_M'] = (df_fit[conc_col] / MW) * 1e-6
        
        t = df_fit['time_s'].values
        E = df_fit['E_M'].values
        y = df_fit['alpha'].values
        
        self._t_data = t
        self._E_data = E
        
        def model_func(x_dummy, kcat_KM, km, Gamma0):
            """Model function for curve_fit"""
            result = np.zeros_like(self._t_data)
            for i in range(len(self._t_data)):
                result[i] = self.model(self._t_data[i], self._E_data[i], kcat_KM, km, Gamma0)
            return result
        
        try:
            # Better initial guess estimation
            kobs_estimates = []
            for conc in df_fit['enzyme_ugml'].unique():
                subset = df_fit[df_fit['enzyme_ugml'] == conc].sort_values('time_s')
                if len(subset) >= 3:
                    early = subset.head(min(5, len(subset)))
                    t_early = early['time_s'].values
                    alpha_early = early['alpha'].values
                    alpha_safe = np.clip(alpha_early, 0.01, 0.99)
                    y_log = -np.log(1 - alpha_safe)
                    if len(t_early) >= 2 and t_early[-1] > t_early[0]:
                        kobs_est = (y_log[-1] - y_log[0]) / (t_early[-1] - t_early[0])
                        E_conc = subset['E_M'].iloc[0]
                        if E_conc > 0 and kobs_est > 0:
                            kcat_KM_est = kobs_est / E_conc
                            kobs_estimates.append(kcat_KM_est)
            
            if len(kobs_estimates) > 0:
                kcat_KM_guess = np.median(kobs_estimates)
                kcat_KM_guess = np.clip(kcat_KM_guess, 1e2, 1e10)
            else:
                kcat_KM_guess = 1e5
            
            p0 = [kcat_KM_guess, 1e-4, 1.0]
            
            if verbose_callback:
                verbose_callback(f"Model C 초기값: kcat/KM = {kcat_KM_guess:.2e} M⁻¹s⁻¹, km = 1e-4 m/s, Γ₀ = 1.0 pmol/cm²")
            
            x_dummy = np.arange(len(y))
            
            popt, pcov = curve_fit(
                model_func, x_dummy, y, p0=p0,
                bounds=([1e3, 1e-6, 0.1], [1e9, 1e-2, 100]),
                maxfev=20000
            )
            
            kcat_KM_fit, km_fit, Gamma0_fit = popt
            perr = np.sqrt(np.diag(pcov))
            
            y_pred = model_func(x_dummy, kcat_KM_fit, km_fit, Gamma0_fit)
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(np.mean(residuals ** 2))
            
            n = len(y)
            k = 3
            if ss_res > 0:
                log_likelihood = -n/2 * np.log(2*np.pi*ss_res/n) - n/2
            else:
                log_likelihood = 0
            aic = 2*k - 2*log_likelihood
            bic = k*np.log(n) - 2*log_likelihood
            
            params = {'kcat_KM': kcat_KM_fit, 'km': km_fit, 'Gamma0': Gamma0_fit}
            params_std = {'kcat_KM': perr[0], 'km': perr[1], 'Gamma0': perr[2]}
            
            t_full = df['time_s'].values
            conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
            
            if 'uM' in conc_col:
                E_full = df[conc_col].values * 1e-6
            elif 'nM' in conc_col:
                E_full = df[conc_col].values * 1e-9
            else:
                E_full = (df[conc_col].values / MW) * 1e-6
            y_pred_full = np.zeros(len(t_full))
            for i in range(len(t_full)):
                y_pred_full[i] = self.model(t_full[i], E_full[i], kcat_KM_fit, km_fit, Gamma0_fit)
            
            return ModelResults(
                name=self.name,
                params=params,
                params_std=params_std,
                aic=aic,
                bic=bic,
                r_squared=r_squared,
                rmse=rmse,
                predictions=y_pred_full,
                residuals=df['alpha'].values - y_pred_full
            )
            
        except Exception as e:
            if verbose_callback:
                verbose_callback(f"Model C 피팅 오류: {str(e)}", level="error")
                import traceback
                verbose_callback(traceback.format_exc(), level="debug")
            return None


class ModelD_ConcentrationDependentFmax:
    """
    Model D: Concentration-Dependent Maximum Cleavage (농도 의존적 Fmax)
    
    Hypothesis: 효소 농도가 높을수록 더 깊은 겔 층까지 침투하여 더 많은 기질을 절단
    
    Physical basis:
    - 확산 침투 깊이: δ ∝ √(D*t) but also depends on reaction rate
    - 높은 [E]에서 더 큰 농도 구배 → 더 깊은 침투
    - 생성물 방출/2차 절단: 높은 [E]에서 작은 조각으로 분해되어 방출 증가
    
    Model equation:
    α(t) = α_max([E]) * [1 - exp(-kobs*t)]
    
    where:
    α_max([E]) = α_∞ * [1 - exp(-k_access * [E])]
    - α_∞: 무한 효소 농도에서 이론적 최대값 (접근 가능한 전체 기질)
    - k_access: 효소 농도에 대한 접근성 계수 (M⁻¹)
    - kobs = (kcat/KM) * [E]
    
    Parameters to fit:
    - kcat_KM: 촉매 효율 (M⁻¹s⁻¹)
    - alpha_inf: 이론적 최대 절단 비율 (0-1)
    - k_access: 농도 의존적 접근성 계수 (M⁻¹)
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        self.enzyme_mw = enzyme_mw
        self.name = "Model D: Concentration-Dependent Fmax"
    
    def model(self, t: np.ndarray, E_M: float, kcat_KM: float, 
              alpha_inf: float, k_access: float) -> np.ndarray:
        """
        Model with concentration-dependent maximum cleavage
        """
        # Concentration-dependent maximum
        alpha_max = alpha_inf * (1.0 - np.exp(-k_access * E_M))
        
        # Time-dependent approach to maximum
        kobs = kcat_KM * E_M
        alpha = alpha_max * (1.0 - np.exp(-kobs * t))
        
        return np.clip(alpha, 0, 1)
    
    def fit_global(self, df: pd.DataFrame, verbose_callback=None) -> ModelResults:
        """Global fit with concentration-dependent Fmax"""
        df_fit = df.copy()
        
        if verbose_callback:
            verbose_callback(f"Model D: 전체 {len(df_fit)}개 데이터 포인트 사용")
        
        if len(df_fit) < 5:
            if verbose_callback:
                verbose_callback(f"Model D: 피팅 데이터가 충분하지 않습니다.", level="error")
            return None
        
        MW = self.enzyme_mw * 1000
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        
        if 'uM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-6
        elif 'nM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-9
        else:
            df_fit['E_M'] = (df_fit[conc_col] / MW) * 1e-6
        
        t = df_fit['time_s'].values
        E = df_fit['E_M'].values
        y = df_fit['alpha'].values
        
        self._t_data = t
        self._E_data = E
        
        def model_func(x_dummy, kcat_KM, alpha_inf, k_access):
            """Model function for curve_fit"""
            result = np.zeros_like(self._t_data)
            for i in range(len(self._t_data)):
                result[i] = self.model(self._t_data[i], self._E_data[i], 
                                      kcat_KM, alpha_inf, k_access)
            return result
        
        try:
            # Initial guess
            kcat_KM_guess = 1e5
            alpha_inf_guess = 1.0  # assume full cleavage is theoretically possible
            k_access_guess = 1e6  # M^-1
            
            p0 = [kcat_KM_guess, alpha_inf_guess, k_access_guess]
            
            if verbose_callback:
                verbose_callback(f"Model D 초기값: kcat/KM={kcat_KM_guess:.2e}, α_∞={alpha_inf_guess:.2f}, k_access={k_access_guess:.2e}")
            
            x_dummy = np.arange(len(y))
            
            popt, pcov = curve_fit(
                model_func, x_dummy, y, p0=p0,
                bounds=([1e3, 0.5, 1e4], [1e9, 1.0, 1e9]),
                maxfev=20000
            )
            
            kcat_KM_fit, alpha_inf_fit, k_access_fit = popt
            perr = np.sqrt(np.diag(pcov))
            
            y_pred = model_func(x_dummy, kcat_KM_fit, alpha_inf_fit, k_access_fit)
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(np.mean(residuals ** 2))
            
            n = len(y)
            k = 3
            if ss_res > 0:
                log_likelihood = -n/2 * np.log(2*np.pi*ss_res/n) - n/2
            else:
                log_likelihood = 0
            aic = 2*k - 2*log_likelihood
            bic = k*np.log(n) - 2*log_likelihood
            
            params = {
                'kcat_KM': kcat_KM_fit,
                'alpha_inf': alpha_inf_fit,
                'k_access': k_access_fit
            }
            params_std = {
                'kcat_KM': perr[0],
                'alpha_inf': perr[1],
                'k_access': perr[2]
            }
            
            t_full = df['time_s'].values
            conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
            
            if 'uM' in conc_col:
                E_full = df[conc_col].values * 1e-6
            elif 'nM' in conc_col:
                E_full = df[conc_col].values * 1e-9
            else:
                E_full = (df[conc_col].values / MW) * 1e-6
            y_pred_full = np.zeros(len(t_full))
            for i in range(len(t_full)):
                y_pred_full[i] = self.model(t_full[i], E_full[i], 
                                           kcat_KM_fit, alpha_inf_fit, k_access_fit)
            
            return ModelResults(
                name=self.name,
                params=params,
                params_std=params_std,
                aic=aic,
                bic=bic,
                r_squared=r_squared,
                rmse=rmse,
                predictions=y_pred_full,
                residuals=df['alpha'].values - y_pred_full
            )
            
        except Exception as e:
            if verbose_callback:
                verbose_callback(f"Model D 피팅 오류: {str(e)}", level="error")
                import traceback
                verbose_callback(traceback.format_exc(), level="debug")
            return None


class ModelE_ProductInhibition:
    """
    Model E: Product Inhibition (생성물 억제)
    
    Hypothesis: 절단된 펩타이드 조각(생성물)이 효소를 경쟁적으로 억제
    
    Physical basis:
    - 생성물이 활성 부위에 결합하여 억제
    - 생성물이 겔/표면에 축적 (느린 확산/방출)
    - 높은 [E]에서는 2차 절단으로 작은 조각 생성 → 억제 감소
    
    Reaction scheme:
    S + E ⇌ ES → P + E
    P + E ⇌ EP (inhibition)
    
    Michaelis-Menten with competitive product inhibition:
    v = (kcat * [E] * [S]) / (KM * (1 + [P]/Ki) + [S])
    
    For [S] << KM (first order):
    dS/dt = -(kcat/KM) * [E] * [S] / (1 + [P]/Ki)
    
    Assume: [P] = [S0] - [S] = [S0] * α
    
    Model equation (simplified):
    dα/dt = kobs * (1 - α) / (1 + Ki_eff * α)
    
    where:
    - kobs = (kcat/KM) * [E]
    - Ki_eff = [S0]/Ki (dimensionless inhibition constant)
    
    Parameters to fit:
    - kcat_KM: catalytic efficiency (M⁻¹s⁻¹)
    - Ki_eff: effective product inhibition constant
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        self.enzyme_mw = enzyme_mw
        self.name = "Model E: Product Inhibition"
    
    def model(self, t: np.ndarray, E_M: float, kcat_KM: float, Ki_eff: float) -> np.ndarray:
        """
        Model with product inhibition
        Numerically integrate: dα/dt = kobs*(1-α)/(1 + Ki_eff*α)
        """
        kobs = kcat_KM * E_M
        
        def dydt(t_val, y):
            alpha = y[0]
            if alpha >= 1.0:
                return [0]
            # Product inhibition term
            rate = kobs * (1 - alpha) / (1 + Ki_eff * alpha)
            return [rate]
        
        sol = solve_ivp(dydt, [0, t[-1] if hasattr(t, '__len__') else t], 
                       [0.0], t_eval=t if hasattr(t, '__len__') else [t], 
                       method='LSODA', rtol=1e-6)
        alpha = sol.y[0]
        return np.clip(alpha, 0, 1)
    
    def fit_global(self, df: pd.DataFrame, verbose_callback=None) -> ModelResults:
        """Global fit with product inhibition"""
        df_fit = df.copy()
        
        if verbose_callback:
            verbose_callback(f"Model E: 전체 {len(df_fit)}개 데이터 포인트 사용")
        
        if len(df_fit) < 5:
            if verbose_callback:
                verbose_callback(f"Model E: 피팅 데이터가 충분하지 않습니다.", level="error")
            return None
        
        MW = self.enzyme_mw * 1000
        conc_col = df_fit['conc_col_name'].iloc[0] if 'conc_col_name' in df_fit.columns else 'enzyme_ugml'
        
        if 'uM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-6
        elif 'nM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-9
        else:
            df_fit['E_M'] = (df_fit[conc_col] / MW) * 1e-6
        
        # Group by concentration for fitting
        conc_list = sorted(df_fit[conc_col].unique())
        
        all_t = []
        all_E = []
        all_y = []
        
        for conc in conc_list:
            subset = df_fit[df_fit[conc_col] == conc].sort_values('time_s')
            all_t.extend(subset['time_s'].values)
            all_E.extend(subset['E_M'].values)
            all_y.extend(subset['alpha'].values)
        
        t = np.array(all_t)
        E = np.array(all_E)
        y = np.array(all_y)
        
        self._conc_list = conc_list
        self._df_fit = df_fit
        self._conc_col = conc_col
        self._MW = MW
        
        def model_func(x_dummy, kcat_KM, Ki_eff):
            """Model function for curve_fit"""
            result = []
            for conc in self._conc_list:
                subset = self._df_fit[self._df_fit[self._conc_col] == conc].sort_values('time_s')
                t_subset = subset['time_s'].values
                E_subset = subset['E_M'].iloc[0]
                
                alpha_pred = self.model(t_subset, E_subset, kcat_KM, Ki_eff)
                result.extend(alpha_pred)
            
            return np.array(result)
        
        try:
            # Initial guess
            kcat_KM_guess = 1e5
            Ki_eff_guess = 1.0
            
            p0 = [kcat_KM_guess, Ki_eff_guess]
            
            if verbose_callback:
                verbose_callback(f"Model E 초기값: kcat/KM={kcat_KM_guess:.2e}, Ki_eff={Ki_eff_guess:.2f}")
            
            x_dummy = np.arange(len(y))
            
            popt, pcov = curve_fit(
                model_func, x_dummy, y, p0=p0,
                bounds=([1e3, 0.01], [1e9, 100]),
                maxfev=20000
            )
            
            kcat_KM_fit, Ki_eff_fit = popt
            perr = np.sqrt(np.diag(pcov))
            
            y_pred = model_func(x_dummy, kcat_KM_fit, Ki_eff_fit)
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(np.mean(residuals ** 2))
            
            n = len(y)
            k = 2
            if ss_res > 0:
                log_likelihood = -n/2 * np.log(2*np.pi*ss_res/n) - n/2
            else:
                log_likelihood = 0
            aic = 2*k - 2*log_likelihood
            bic = k*np.log(n) - 2*log_likelihood
            
            params = {'kcat_KM': kcat_KM_fit, 'Ki_eff': Ki_eff_fit}
            params_std = {'kcat_KM': perr[0], 'Ki_eff': perr[1]}
            
            # Generate predictions for all data
            y_pred_full = []
            conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
            
            for conc in sorted(df[conc_col].unique()):
                subset = df[df[conc_col] == conc].sort_values('time_s')
                t_subset = subset['time_s'].values
                
                # Convert concentration to M
                if 'uM' in conc_col:
                    E_subset = conc * 1e-6
                elif 'nM' in conc_col:
                    E_subset = conc * 1e-9
                else:
                    E_subset = (conc / MW) * 1e-6
                
                alpha_pred = self.model(t_subset, E_subset, kcat_KM_fit, Ki_eff_fit)
                y_pred_full.extend(alpha_pred)
            
            y_pred_full = np.array(y_pred_full)
            
            return ModelResults(
                name=self.name,
                params=params,
                params_std=params_std,
                aic=aic,
                bic=bic,
                r_squared=r_squared,
                rmse=rmse,
                predictions=y_pred_full,
                residuals=df['alpha'].values - y_pred_full
            )
            
        except Exception as e:
            if verbose_callback:
                verbose_callback(f"Model E 피팅 오류: {str(e)}", level="error")
                import traceback
                verbose_callback(traceback.format_exc(), level="debug")
            return None


class ModelF_EnzymeSurfaceSequestration:
    """
    Model F: Enzyme Surface Sequestration/Adsorption (효소 표면 흡착/격리)
    
    Hypothesis: 효소가 겔/표면에 비특이적으로 흡착되어 비가역적으로 비활성화
    
    Physical basis:
    - PDA·음전하 표면에 양전하 효소 흡착
    - 겔 망상구조에 물리적 포획
    - 높은 [E]에서 포화/경쟁으로 상대적 영향 감소
    
    Model:
    E_free(t) = E0 * exp(-k_ads * t) / (1 + K_ads * E0)
    
    Or simplified two-state:
    - Fast reversible binding: E ⇌ E_bound (K_eq)
    - Slow irreversible adsorption: E → E_ads (k_ads)
    
    Effective enzyme:
    E_eff([E], t) = [E] * exp(-k_ads * t) / (1 + K_ads * [E])
    
    Integrated:
    α(t) = ∫[0,t] (kcat/KM) * E_eff(τ) dτ
    
    Approximate solution:
    α(t) ≈ (kcat/KM) * [E] / (k_ads * (1 + K_ads*[E])) * [1 - exp(-k_ads*t)]
    
    Parameters to fit:
    - kcat_KM: catalytic efficiency (M⁻¹s⁻¹)
    - k_ads: adsorption rate constant (s⁻¹)
    - K_ads: adsorption equilibrium constant (M⁻¹)
    """
    
    def __init__(self, enzyme_mw: float = 56.6):
        self.enzyme_mw = enzyme_mw
        self.name = "Model F: Enzyme Surface Sequestration"
    
    def model(self, t: np.ndarray, E_M: float, kcat_KM: float, 
              k_ads: float, K_ads: float) -> np.ndarray:
        """
        Model with enzyme surface adsorption/sequestration
        """
        # Concentration-dependent adsorption factor
        ads_factor = 1.0 / (1.0 + K_ads * E_M)
        
        # Time-dependent depletion with saturation
        if k_ads > 1e-6:
            alpha = (kcat_KM * E_M * ads_factor / k_ads) * (1.0 - np.exp(-k_ads * t))
        else:
            # Limit as k_ads → 0 (no adsorption)
            alpha = kcat_KM * E_M * ads_factor * t
        
        return np.clip(alpha, 0, 1)
    
    def fit_global(self, df: pd.DataFrame, verbose_callback=None) -> ModelResults:
        """Global fit with enzyme sequestration"""
        df_fit = df.copy()
        
        if verbose_callback:
            verbose_callback(f"Model F: 전체 {len(df_fit)}개 데이터 포인트 사용")
        
        if len(df_fit) < 5:
            if verbose_callback:
                verbose_callback(f"Model F: 피팅 데이터가 충분하지 않습니다.", level="error")
            return None
        
        MW = self.enzyme_mw * 1000
        conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
        
        if 'uM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-6
        elif 'nM' in conc_col:
            df_fit['E_M'] = df_fit[conc_col] * 1e-9
        else:
            df_fit['E_M'] = (df_fit[conc_col] / MW) * 1e-6
        
        t = df_fit['time_s'].values
        E = df_fit['E_M'].values
        y = df_fit['alpha'].values
        
        self._t_data = t
        self._E_data = E
        
        def model_func(x_dummy, kcat_KM, k_ads, K_ads):
            """Model function for curve_fit"""
            result = np.zeros_like(self._t_data)
            for i in range(len(self._t_data)):
                result[i] = self.model(self._t_data[i], self._E_data[i], 
                                      kcat_KM, k_ads, K_ads)
            return result
        
        try:
            # Initial guess
            kcat_KM_guess = 1e5
            k_ads_guess = 0.1  # s^-1
            K_ads_guess = 1e6  # M^-1
            
            p0 = [kcat_KM_guess, k_ads_guess, K_ads_guess]
            
            if verbose_callback:
                verbose_callback(f"Model F 초기값: kcat/KM={kcat_KM_guess:.2e}, k_ads={k_ads_guess:.2f}, K_ads={K_ads_guess:.2e}")
            
            x_dummy = np.arange(len(y))
            
            popt, pcov = curve_fit(
                model_func, x_dummy, y, p0=p0,
                bounds=([1e3, 0.001, 1e3], [1e9, 10, 1e9]),
                maxfev=20000
            )
            
            kcat_KM_fit, k_ads_fit, K_ads_fit = popt
            perr = np.sqrt(np.diag(pcov))
            
            y_pred = model_func(x_dummy, kcat_KM_fit, k_ads_fit, K_ads_fit)
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(np.mean(residuals ** 2))
            
            n = len(y)
            k = 3
            if ss_res > 0:
                log_likelihood = -n/2 * np.log(2*np.pi*ss_res/n) - n/2
            else:
                log_likelihood = 0
            aic = 2*k - 2*log_likelihood
            bic = k*np.log(n) - 2*log_likelihood
            
            params = {
                'kcat_KM': kcat_KM_fit,
                'k_ads': k_ads_fit,
                'K_ads': K_ads_fit
            }
            params_std = {
                'kcat_KM': perr[0],
                'k_ads': perr[1],
                'K_ads': perr[2]
            }
            
            t_full = df['time_s'].values
            conc_col = df['conc_col_name'].iloc[0] if 'conc_col_name' in df.columns else 'enzyme_ugml'
            
            if 'uM' in conc_col:
                E_full = df[conc_col].values * 1e-6
            elif 'nM' in conc_col:
                E_full = df[conc_col].values * 1e-9
            else:
                E_full = (df[conc_col].values / MW) * 1e-6
            y_pred_full = np.zeros(len(t_full))
            for i in range(len(t_full)):
                y_pred_full[i] = self.model(t_full[i], E_full[i], 
                                           kcat_KM_fit, k_ads_fit, K_ads_fit)
            
            return ModelResults(
                name=self.name,
                params=params,
                params_std=params_std,
                aic=aic,
                bic=bic,
                r_squared=r_squared,
                rmse=rmse,
                predictions=y_pred_full,
                residuals=df['alpha'].values - y_pred_full
            )
            
        except Exception as e:
            if verbose_callback:
                verbose_callback(f"Model F 피팅 오류: {str(e)}", level="error")
                import traceback
                verbose_callback(traceback.format_exc(), level="debug")
            return None

