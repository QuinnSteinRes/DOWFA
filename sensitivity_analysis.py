"""
Sensitivity Analysis for Volatile Anaesthetic Waste Disposal Routes
--------------------------------------------------------------------
Implements Standardized Regression Coefficients (SRC) and 
Partial Rank Correlation Coefficients (PRCC) for routes R3, R4, R5.

Based on reviewer response R2-3 requirements.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from parameters_m import *

# ============================================================================
# SENSITIVITY ANALYSIS FUNCTIONS
# ============================================================================

def _standardize(a, axis=0):
    """Z-score standardization (mean=0, std=1)"""
    a = np.asarray(a, dtype=float)
    m = a.mean(axis=axis, keepdims=True)
    s = a.std(axis=axis, ddof=1, keepdims=True)
    return (a - m) / s, m.squeeze(), s.squeeze()

def _ols_beta(X, y):
    """Ordinary least squares coefficients"""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta

def src_coefficients(X, y):
    """
    Standardized Regression Coefficients (SRC)
    
    Fits OLS on z-scored inputs/outputs:
    ỹ = α + Σ βⱼ(S) X̃ⱼ + ε
    
    Returns the standardized coefficients βⱼ(S)
    """
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    Xs, _, _ = _standardize(X, axis=0)
    ys, _, _ = _standardize(y, axis=0)
    Z = np.column_stack([np.ones(Xs.shape[0]), Xs])
    b = _ols_beta(Z, ys)
    return b[1:]  # drop intercept

def prcc_coefficients(X, y, rank_method="average"):
    """
    Partial Rank Correlation Coefficients (PRCC)
    
    For each input Xj:
    1. Rank all variables
    2. Regress rank(Xj) on rank(other inputs)
    3. Regress rank(y) on rank(other inputs)  
    4. Compute correlation of residuals
    """
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    n, p = X.shape
    
    def rank_col(v):
        return pd.Series(v).rank(method=rank_method).to_numpy()
    
    XR = np.column_stack([rank_col(X[:, j]) for j in range(p)])
    yR = rank_col(y)
    prcc = np.zeros(p, dtype=float)
    
    for j in range(p):
        others = [i for i in range(p) if i != j]
        W = np.column_stack([np.ones(n), XR[:, others]])
        
        # Residuals: Xj | X_-j
        bj = _ols_beta(W, XR[:, j])
        e_j = XR[:, j] - W @ bj
        
        # Residuals: y | X_-j  
        by = _ols_beta(W, yR)
        e_y = yR - W @ by
        
        prcc[j] = np.corrcoef(e_j, e_y)[0, 1]
    
    return prcc

def tornado_plot(names, values, title, xlabel, outfile=None):
    """
    Horizontal bar chart ordered by |values|
    """
    order = np.argsort(np.abs(values))
    names_o = [names[i] for i in order]
    vals_o = [values[i] for i in order]
    
    plt.figure(figsize=(8, 5))
    y = np.arange(len(vals_o))
    plt.barh(y, vals_o, color='steelblue')
    plt.yticks(y, names_o)
    plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.show()

# ============================================================================
# ROUTE MODELS (per kg VA)
# ============================================================================

def R3_MSWI(x3, x4, x5, x6, anaesthetic='des', gwp_type='100'):
    """
    R3: Municipal Solid Waste Incineration (MSWI)
    
    Parameters from Table 2:
    x3: HGV factor [kgCO2e/tonne/km]
    x4: Distance to MSWI [km]
    x5: CF4 generation fraction [0..1]
    x6: Generic combustion emissions [kgCO2e/kgVA]
    
    Model: R3 = (x3·x4)/1000 + (1/MW_VA)[formation terms] + x6
    """
    # Set constants based on anaesthetic and GWP
    if anaesthetic == 'des':
        MW_VA = des_gfm / 1000  # Convert to kg/mol
        carbon_atoms = des_atecon_C
    elif anaesthetic == 'sev':
        MW_VA = sev_gfm / 1000
        carbon_atoms = sev_atecon_C
    elif anaesthetic == 'iso':
        MW_VA = iso_gfm / 1000
        carbon_atoms = iso_atecon_C
    
    # GWP values
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044  # kg/mol
    MW_CF4 = 0.088  # kg/mol
    
    # Transport term
    transport = (x3 * x4) / 1000.0
    
    # Formation terms
    term_CO2 = ((1.0 - x5) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (x5 * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    return transport + term_CO2 + term_CF4 + x6

def R4_HCWI(x3, x4, x7, x8, anaesthetic='des', gwp_type='100'):
    """
    R4: Healthcare Waste Incineration (HCWI)
    
    Parameters from Table 2:
    x3: HGV factor [kgCO2e/tonne/km]
    x4: Distance to HCWI [km] (using same as MSWI for simplicity)
    x7: CF4 generation fraction in HCWI [0..1]
    x8: Generic combustion emissions HCWI [kgCO2e/kgVA]
    
    Model: Similar to R3 but with HCWI-specific parameters
    """
    # Set constants
    if anaesthetic == 'des':
        MW_VA = des_gfm / 1000
        carbon_atoms = des_atecon_C
    elif anaesthetic == 'sev':
        MW_VA = sev_gfm / 1000
        carbon_atoms = sev_atecon_C
    elif anaesthetic == 'iso':
        MW_VA = iso_gfm / 1000
        carbon_atoms = iso_atecon_C
    
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044
    MW_CF4 = 0.088
    
    transport = (x3 * x4) / 1000.0
    term_CO2 = ((1.0 - x7) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (x7 * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    return transport + term_CO2 + term_CF4 + x8

def R5_plasma(x9, anaesthetic='des', gwp_type='100'):
    """
    R5: In-situ Plasma Destruction
    
    Parameters from Table 2:
    x9: Plasma destruction emissions [kgCO2e/kgVA]
    
    Model: R5 = x9 (very simple - just the destruction technology emissions)
    """
    return x9

# ============================================================================
# SENSITIVITY ANALYSIS RUNNERS
# ============================================================================

def run_R3_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R3 (MSWI)
    
    Varied inputs (from Table 2):
    - x3: HGV factor [0.01816, 0.12465] kgCO2e/tonne/km
    - x4: Distance to MSWI [24, 48] km
    - x5: CF4 generation fraction [0.17, 0.40]
    - x6: Combustion CO2e [0.7, 1.7] kgCO2e/kgVA
    """
    rng = np.random.default_rng(seed)
    
    # Sample from uniform ranges
    x3 = rng.uniform(min_hgv_factor, max_hgv_factor, N)
    x4 = rng.uniform(min_mwi_dist, max_mwi_dist, N)
    x5 = rng.uniform(min_demunic, max_demunic, N)
    x6 = rng.uniform(min_munic_other, max_munic_other, N)
    
    # Evaluate model
    y = R3_MSWI(x3, x4, x5, x6, anaesthetic, gwp_type)
    
    # Build input matrix
    X = np.column_stack([x3, x4, x5, x6])
    names = ["HGV factor (x3)", "Distance to MSWI (x4)", 
             "CF4 fraction (x5)", "Combustion CO2e (x6)"]
    
    # Compute sensitivity metrics
    SRC = src_coefficients(X, y)
    PRCC = prcc_coefficients(X, y)
    
    # Create results table
    df = pd.DataFrame({
        "Input": names,
        "SRC": SRC,
        "PRCC": PRCC,
        "Abs(SRC)": np.abs(SRC),
        "Abs(PRCC)": np.abs(PRCC),
    }).sort_values("Abs(SRC)", ascending=False).reset_index(drop=True)
    
    # Summary statistics
    summary = {
        "route": "R3_MSWI",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    # Save results
    if save_csv:
        filename = f"{outprefix}R3_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    # Create tornado plot
    if make_plot:
        plot_file = f"{outprefix}R3_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R3 MSWI Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_R4_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R4 (HCWI)
    
    Varied inputs:
    - x3: HGV factor [0.01816, 0.12465] kgCO2e/tonne/km
    - x4: Distance [24, 48] km (using MSWI distance range)
    - x7: CF4 generation fraction [0.17, 0.40] (assuming similar to MSWI)
    - x8: Combustion CO2e [0.851, 1.05] kgCO2e/kgVA
    """
    rng = np.random.default_rng(seed)
    
    x3 = rng.uniform(min_hgv_factor, max_hgv_factor, N)
    x4 = rng.uniform(min_mwi_dist, max_mwi_dist, N)
    x7 = rng.uniform(min_demunic, max_demunic, N)  # Using MSWI range
    x8 = rng.uniform(min_pharam_other, max_pharam_other, N)
    
    y = R4_HCWI(x3, x4, x7, x8, anaesthetic, gwp_type)
    
    X = np.column_stack([x3, x4, x7, x8])
    names = ["HGV factor (x3)", "Distance to HCWI (x4)",
             "CF4 fraction (x7)", "Combustion CO2e (x8)"]
    
    SRC = src_coefficients(X, y)
    PRCC = prcc_coefficients(X, y)
    
    df = pd.DataFrame({
        "Input": names,
        "SRC": SRC,
        "PRCC": PRCC,
        "Abs(SRC)": np.abs(SRC),
        "Abs(PRCC)": np.abs(PRCC),
    }).sort_values("Abs(SRC)", ascending=False).reset_index(drop=True)
    
    summary = {
        "route": "R4_HCWI",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R4_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    if make_plot:
        plot_file = f"{outprefix}R4_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R4 HCWI Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_R5_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R5 (Plasma destruction)
    
    Varied input:
    - x9: Plasma destruction emissions [0.040, 0.089] kgCO2e/kgVA
    
    Note: R5 has only one parameter, so sensitivity analysis is trivial
    (perfect correlation), but included for completeness.
    """
    rng = np.random.default_rng(seed)
    
    x9 = rng.uniform(min_plasma_other, max_plasma_other, N)
    y = R5_plasma(x9, anaesthetic, gwp_type)
    
    X = np.column_stack([x9])
    names = ["Plasma emissions (x9)"]
    
    SRC = src_coefficients(X, y)
    PRCC = prcc_coefficients(X, y)
    
    df = pd.DataFrame({
        "Input": names,
        "SRC": SRC,
        "PRCC": PRCC,
        "Abs(SRC)": np.abs(SRC),
        "Abs(PRCC)": np.abs(PRCC),
    })
    
    summary = {
        "route": "R5_Plasma",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R5_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    print("\nNote: R5 has only one parameter, so SRC/PRCC are trivially ±1.0")
    
    return df, summary

def run_all_routes_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                               save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for all three routes (R3, R4, R5)
    """
    print("="*80)
    print(f"SENSITIVITY ANALYSIS: {anaesthetic.upper()}, GWP-{gwp_type}")
    print("="*80)
    
    print("\n--- Route R3: MSWI ---")
    df_r3, sum_r3 = run_R3_sensitivity(seed, N, anaesthetic, gwp_type, 
                                       save_csv, make_plot, outprefix)
    print(df_r3)
    print(f"\nOutput mean: {sum_r3['mean']:.2f} ± {sum_r3['std']:.2f} kgCO2e/kg")
    
    print("\n--- Route R4: HCWI ---")
    df_r4, sum_r4 = run_R4_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r4)
    print(f"\nOutput mean: {sum_r4['mean']:.2f} ± {sum_r4['std']:.2f} kgCO2e/kg")
    
    print("\n--- Route R5: Plasma ---")
    df_r5, sum_r5 = run_R5_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r5)
    print(f"\nOutput mean: {sum_r5['mean']:.2f} ± {sum_r5['std']:.2f} kgCO2e/kg")
    
    print("\n" + "="*80)
    
    return {
        'R3': (df_r3, sum_r3),
        'R4': (df_r4, sum_r4),
        'R5': (df_r5, sum_r5)
    }

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Run for desflurane, GWP-100 (as in the LaTeX example)
    results = run_all_routes_sensitivity(
        seed=42, 
        N=10000,
        anaesthetic='des',
        gwp_type='100',
        save_csv=True,
        make_plot=True,
        outprefix="./"
    )
