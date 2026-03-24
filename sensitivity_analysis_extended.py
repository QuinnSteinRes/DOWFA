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

# Routes 6-11: Waste Volatile Anaesthetic Routes
# These involve activated carbon (AC) capture followed by various disposal methods

def R6_continuous_release(anaesthetic='des', gwp_type='100'):
    """
    R6: Continuous atmospheric release (baseline/reference case)
    
    This is just the direct GWP of the anaesthetic released to atmosphere.
    No parameters vary - this is deterministic given the anaesthetic and GWP choice.
    
    Returns the GWP value for 1 kg of the anaesthetic.
    """
    if gwp_type == '20':
        gwp_map = {'des': des_20, 'sev': sev_20, 'iso': iso_20}
    elif gwp_type == '100':
        gwp_map = {'des': des_100, 'sev': sev_100, 'iso': iso_100}
    elif gwp_type == '500':
        gwp_map = {'des': des_500, 'sev': sev_500, 'iso': iso_500}
    
    return np.ones(1) * gwp_map[anaesthetic]  # Return as array for consistency

def R7_AC_MSWI(van_f, hgv_f, ship_f, net_dist, mwi_dist, cf4_frac, combust_co2e, 
               anaesthetic='des', gwp_type='100'):
    """
    R7: AC capture + transport to MSWI + incineration
    
    Parameters:
    van_f: Van transport factor [kgCO2e/tonne/km]
    hgv_f: HGV transport factor [kgCO2e/tonne/km]
    ship_f: Ship transport factor [kgCO2e/tonne/km]
    net_dist: Network distribution distance [km]
    mwi_dist: Distance to MSWI [km]
    cf4_frac: CF4 generation fraction at MSWI [0..1]
    combust_co2e: Generic combustion emissions [kgCO2e/kgVA]
    
    Model components:
    1. AC production emissions (per kg AC)
    2. AC transport (Indonesia→Singapore→UK→facility)
    3. Network distribution (facility→hospitals)
    4. Final transport to MSWI (with loaded AC + VA)
    5. Combustion emissions (CO2 + CF4 formation)
    6. Generic incineration emissions
    """
    # Get anaesthetic properties
    if anaesthetic == 'des':
        mol_per_kg_ac = mol_per_kg_des
        carbon_atoms = des_atecon_C
        MW_VA = des_gfm / 1000
    elif anaesthetic == 'sev':
        mol_per_kg_ac = mol_per_kg_sev
        carbon_atoms = sev_atecon_C
        MW_VA = sev_gfm / 1000
    elif anaesthetic == 'iso':
        mol_per_kg_ac = mol_per_kg_iso
        carbon_atoms = iso_atecon_C
        MW_VA = iso_gfm / 1000
    
    # GWP
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044
    MW_CF4 = 0.088
    
    # AC requirement: kg AC needed per kg VA
    kg_ac_per_kg_va = 1.0 / (mol_per_kg_ac * MW_VA)
    
    # 1. AC production: assume coconut ratio and production emission factor
    ac_prod = kg_ac_per_kg_va * ac_co2e_per_kg
    
    # 2. AC transport (international + to facility)
    # Simplified: Indo→Sing→UK supply chain
    ac_transport_intl = kg_ac_per_kg_va * (
        (indo_to_sing_port * cococnut_to_ac_ratio * ship_f / 1000) +  # Coconut
        (sing_port_to_fac * cococnut_to_ac_ratio * hgv_f / 1000) +     # Coconut to AC facility
        (sing_port_to_fac * hgv_f / 1000) +                            # AC locally
        (sing_to_eng_port * ship_f / 1000)                             # AC to UK
    )
    
    # 3. Network distribution
    ac_network = kg_ac_per_kg_va * net_dist * hgv_f / 1000
    
    # 4. Transport to MSWI (AC + VA together)
    transport_to_mswi = (kg_ac_per_kg_va + 1.0) * mwi_dist * hgv_f / 1000
    
    # 5. Combustion (CO2 and CF4 formation)
    term_CO2 = ((1.0 - cf4_frac) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (cf4_frac * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    # 6. Generic combustion/incineration emissions
    # Applied to both AC and VA mass
    generic_combust = (kg_ac_per_kg_va + 1.0) * combust_co2e
    
    total = ac_prod + ac_transport_intl + ac_network + transport_to_mswi + term_CO2 + term_CF4 + generic_combust
    return total

def R8_AC_HCWI(van_f, hgv_f, ship_f, net_dist, pwi_dist, cf4_frac, combust_co2e,
               anaesthetic='des', gwp_type='100'):
    """
    R8: AC capture + transport to HCWI (pharmaceutical waste incinerator)
    
    Same structure as R7 but with HCWI-specific parameters.
    """
    # Get anaesthetic properties
    if anaesthetic == 'des':
        mol_per_kg_ac = mol_per_kg_des
        carbon_atoms = des_atecon_C
        MW_VA = des_gfm / 1000
    elif anaesthetic == 'sev':
        mol_per_kg_ac = mol_per_kg_sev
        carbon_atoms = sev_atecon_C
        MW_VA = sev_gfm / 1000
    elif anaesthetic == 'iso':
        mol_per_kg_ac = mol_per_kg_iso
        carbon_atoms = iso_atecon_C
        MW_VA = iso_gfm / 1000
    
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044
    MW_CF4 = 0.088
    
    kg_ac_per_kg_va = 1.0 / (mol_per_kg_ac * MW_VA)
    
    ac_prod = kg_ac_per_kg_va * ac_co2e_per_kg
    
    ac_transport_intl = kg_ac_per_kg_va * (
        (indo_to_sing_port * cococnut_to_ac_ratio * ship_f / 1000) +
        (sing_port_to_fac * cococnut_to_ac_ratio * hgv_f / 1000) +
        (sing_port_to_fac * hgv_f / 1000) +
        (sing_to_eng_port * ship_f / 1000)
    )
    
    ac_network = kg_ac_per_kg_va * net_dist * hgv_f / 1000
    
    transport_to_hcwi = (kg_ac_per_kg_va + 1.0) * pwi_dist * hgv_f / 1000
    
    term_CO2 = ((1.0 - cf4_frac) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (cf4_frac * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    generic_combust = (kg_ac_per_kg_va + 1.0) * combust_co2e
    
    total = ac_prod + ac_transport_intl + ac_network + transport_to_hcwi + term_CO2 + term_CF4 + generic_combust
    return total

def R9_AC_plasma(van_f, hgv_f, ship_f, net_dist, cf4_frac, plasma_co2e,
                 anaesthetic='des', gwp_type='100'):
    """
    R9: AC capture + transport for plasma destruction
    
    Similar to R7/R8 but final stage is plasma destruction instead of incineration.
    """
    if anaesthetic == 'des':
        mol_per_kg_ac = mol_per_kg_des
        carbon_atoms = des_atecon_C
        MW_VA = des_gfm / 1000
    elif anaesthetic == 'sev':
        mol_per_kg_ac = mol_per_kg_sev
        carbon_atoms = sev_atecon_C
        MW_VA = sev_gfm / 1000
    elif anaesthetic == 'iso':
        mol_per_kg_ac = mol_per_kg_iso
        carbon_atoms = iso_atecon_C
        MW_VA = iso_gfm / 1000
    
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044
    MW_CF4 = 0.088
    
    kg_ac_per_kg_va = 1.0 / (mol_per_kg_ac * MW_VA)
    
    ac_prod = kg_ac_per_kg_va * ac_co2e_per_kg
    
    ac_transport_intl = kg_ac_per_kg_va * (
        (indo_to_sing_port * cococnut_to_ac_ratio * ship_f / 1000) +
        (sing_port_to_fac * cococnut_to_ac_ratio * hgv_f / 1000) +
        (sing_port_to_fac * hgv_f / 1000) +
        (sing_to_eng_port * ship_f / 1000)
    )
    
    ac_network = kg_ac_per_kg_va * net_dist * hgv_f / 1000
    
    # Transport for plasma (network distance)
    transport_to_plasma = (kg_ac_per_kg_va + 1.0) * net_dist * hgv_f / 1000
    
    # Plasma destruction products
    term_CO2 = ((1.0 - cf4_frac) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (cf4_frac * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    # Plasma destruction emissions (electrical energy)
    plasma_emissions = (kg_ac_per_kg_va + 1.0) * plasma_co2e
    
    total = ac_prod + ac_transport_intl + ac_network + transport_to_plasma + term_CO2 + term_CF4 + plasma_emissions
    return total

def R10_AC_onsite_HCWI(van_f, hgv_f, ship_f, net_dist, cf4_frac, combust_co2e,
                       anaesthetic='des', gwp_type='100'):
    """
    R10: AC capture + on-site incineration at hospital
    
    No transport to external facility, but still need AC production/delivery.
    """
    if anaesthetic == 'des':
        mol_per_kg_ac = mol_per_kg_des
        carbon_atoms = des_atecon_C
        MW_VA = des_gfm / 1000
    elif anaesthetic == 'sev':
        mol_per_kg_ac = mol_per_kg_sev
        carbon_atoms = sev_atecon_C
        MW_VA = sev_gfm / 1000
    elif anaesthetic == 'iso':
        mol_per_kg_ac = mol_per_kg_iso
        carbon_atoms = iso_atecon_C
        MW_VA = iso_gfm / 1000
    
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044
    MW_CF4 = 0.088
    
    kg_ac_per_kg_va = 1.0 / (mol_per_kg_ac * MW_VA)
    
    ac_prod = kg_ac_per_kg_va * ac_co2e_per_kg
    
    ac_transport_intl = kg_ac_per_kg_va * (
        (indo_to_sing_port * cococnut_to_ac_ratio * ship_f / 1000) +
        (sing_port_to_fac * cococnut_to_ac_ratio * hgv_f / 1000) +
        (sing_port_to_fac * hgv_f / 1000) +
        (sing_to_eng_port * ship_f / 1000)
    )
    
    ac_network = kg_ac_per_kg_va * net_dist * hgv_f / 1000
    
    # No external transport - on-site
    
    # Combustion
    term_CO2 = ((1.0 - cf4_frac) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (cf4_frac * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    generic_combust = (kg_ac_per_kg_va + 1.0) * combust_co2e
    
    total = ac_prod + ac_transport_intl + ac_network + term_CO2 + term_CF4 + generic_combust
    return total

def R11_insitu_plasma(cf4_frac, plasma_co2e, anaesthetic='des', gwp_type='100'):
    """
    R11: In-situ plasma destruction (no AC capture)
    
    Direct destruction at point of generation using plasma technology.
    """
    if anaesthetic == 'des':
        carbon_atoms = des_atecon_C
        MW_VA = des_gfm / 1000
    elif anaesthetic == 'sev':
        carbon_atoms = sev_atecon_C
        MW_VA = sev_gfm / 1000
    elif anaesthetic == 'iso':
        carbon_atoms = iso_atecon_C
        MW_VA = iso_gfm / 1000
    
    if gwp_type == '20':
        GWP_CF4 = cf4_20
    elif gwp_type == '100':
        GWP_CF4 = cf4_100
    elif gwp_type == '500':
        GWP_CF4 = cf4_500
    
    GWP_CO2 = 1.0
    MW_CO2 = 0.044
    MW_CF4 = 0.088
    
    # Destruction products
    term_CO2 = ((1.0 - cf4_frac) * carbon_atoms * GWP_CO2 * MW_CO2) / MW_VA
    term_CF4 = (cf4_frac * carbon_atoms * GWP_CF4 * MW_CF4) / MW_VA
    
    # Plasma emissions
    total = term_CO2 + term_CF4 + plasma_co2e
    return total

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

def run_R6_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=False, outprefix="./"):
    """
    Run sensitivity analysis for R6 (Continuous release)
    
    R6 has NO varying parameters - it's just the GWP of the anaesthetic.
    This is a reference/baseline case. No sensitivity analysis needed.
    """
    y = R6_continuous_release(anaesthetic, gwp_type)
    
    summary = {
        "route": "R6_Continuous_Release",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y[0]),
        "std": 0.0,
        "median": float(y[0]),
        "N": 1,
    }
    
    print(f"\nR6 is deterministic: {y[0]:.1f} kgCO2e/kg (no uncertainty)")
    
    return None, summary

def run_R7_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R7 (AC capture + MSWI)
    
    Varied inputs:
    - van_factor: [0.13926, 0.29416] kgCO2e/tonne/km
    - hgv_factor: [0.01816, 0.12465] kgCO2e/tonne/km
    - ship_factor: [0.00057, 0.00833] kgCO2e/tonne/km
    - network_distance: [2500, 3500] km
    - mwi_dist: [24, 48] km
    - cf4_frac: [0.17, 0.40]
    - combust_co2e: [0.7, 1.7] kgCO2e/kgVA
    """
    rng = np.random.default_rng(seed)
    
    van_f = rng.uniform(min_van_factor, max_van_factor, N)
    hgv_f = rng.uniform(min_hgv_factor, max_hgv_factor, N)
    ship_f = rng.uniform(min_ship_factor, max_ship_factor, N)
    net_dist = rng.uniform(min_network_distance, max_network_distance, N)
    mwi_dist = rng.uniform(min_mwi_dist, max_mwi_dist, N)
    cf4_frac = rng.uniform(min_demunic, max_demunic, N)
    combust_co2e = rng.uniform(min_munic_other, max_munic_other, N)
    
    y = R7_AC_MSWI(van_f, hgv_f, ship_f, net_dist, mwi_dist, cf4_frac, combust_co2e,
                   anaesthetic, gwp_type)
    
    X = np.column_stack([van_f, hgv_f, ship_f, net_dist, mwi_dist, cf4_frac, combust_co2e])
    names = ["Van factor", "HGV factor", "Ship factor", "Network distance", 
             "MSWI distance", "CF4 fraction", "Combustion CO2e"]
    
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
        "route": "R7_AC_MSWI",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R7_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    if make_plot:
        plot_file = f"{outprefix}R7_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R7 AC+MSWI Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_R8_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R8 (AC capture + HCWI)
    """
    rng = np.random.default_rng(seed)
    
    van_f = rng.uniform(min_van_factor, max_van_factor, N)
    hgv_f = rng.uniform(min_hgv_factor, max_hgv_factor, N)
    ship_f = rng.uniform(min_ship_factor, max_ship_factor, N)
    net_dist = rng.uniform(min_network_distance, max_network_distance, N)
    pwi_dist = rng.uniform(min_pwi_dist, max_pwi_dist, N)
    cf4_frac = rng.uniform(min_depharma, max_depharma, N)
    combust_co2e = rng.uniform(min_pharam_other, max_pharam_other, N)
    
    y = R8_AC_HCWI(van_f, hgv_f, ship_f, net_dist, pwi_dist, cf4_frac, combust_co2e,
                   anaesthetic, gwp_type)
    
    X = np.column_stack([van_f, hgv_f, ship_f, net_dist, pwi_dist, cf4_frac, combust_co2e])
    names = ["Van factor", "HGV factor", "Ship factor", "Network distance",
             "HCWI distance", "CF4 fraction", "Combustion CO2e"]
    
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
        "route": "R8_AC_HCWI",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R8_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    if make_plot:
        plot_file = f"{outprefix}R8_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R8 AC+HCWI Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_R9_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                      save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R9 (AC capture + plasma destruction)
    """
    rng = np.random.default_rng(seed)
    
    van_f = rng.uniform(min_van_factor, max_van_factor, N)
    hgv_f = rng.uniform(min_hgv_factor, max_hgv_factor, N)
    ship_f = rng.uniform(min_ship_factor, max_ship_factor, N)
    net_dist = rng.uniform(min_network_distance, max_network_distance, N)
    cf4_frac = rng.uniform(min_deplasma, max_deplasma, N)
    plasma_co2e = rng.uniform(min_plasma_other, max_plasma_other, N)
    
    y = R9_AC_plasma(van_f, hgv_f, ship_f, net_dist, cf4_frac, plasma_co2e,
                     anaesthetic, gwp_type)
    
    X = np.column_stack([van_f, hgv_f, ship_f, net_dist, cf4_frac, plasma_co2e])
    names = ["Van factor", "HGV factor", "Ship factor", "Network distance",
             "CF4 fraction", "Plasma CO2e"]
    
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
        "route": "R9_AC_Plasma",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R9_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    if make_plot:
        plot_file = f"{outprefix}R9_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R9 AC+Plasma Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_R10_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                       save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R10 (AC capture + on-site incineration)
    """
    rng = np.random.default_rng(seed)
    
    van_f = rng.uniform(min_van_factor, max_van_factor, N)
    hgv_f = rng.uniform(min_hgv_factor, max_hgv_factor, N)
    ship_f = rng.uniform(min_ship_factor, max_ship_factor, N)
    net_dist = rng.uniform(min_network_distance, max_network_distance, N)
    cf4_frac = rng.uniform(min_depharma, max_depharma, N)
    combust_co2e = rng.uniform(min_pharam_other, max_pharam_other, N)
    
    y = R10_AC_onsite_HCWI(van_f, hgv_f, ship_f, net_dist, cf4_frac, combust_co2e,
                           anaesthetic, gwp_type)
    
    X = np.column_stack([van_f, hgv_f, ship_f, net_dist, cf4_frac, combust_co2e])
    names = ["Van factor", "HGV factor", "Ship factor", "Network distance",
             "CF4 fraction", "Combustion CO2e"]
    
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
        "route": "R10_AC_Onsite",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R10_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    if make_plot:
        plot_file = f"{outprefix}R10_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R10 AC+Onsite Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_R11_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                       save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for R11 (In-situ plasma destruction)
    """
    rng = np.random.default_rng(seed)
    
    cf4_frac = rng.uniform(min_deplasma, max_deplasma, N)
    plasma_co2e = rng.uniform(min_plasma_other, max_plasma_other, N)
    
    y = R11_insitu_plasma(cf4_frac, plasma_co2e, anaesthetic, gwp_type)
    
    X = np.column_stack([cf4_frac, plasma_co2e])
    names = ["CF4 fraction", "Plasma CO2e"]
    
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
        "route": "R11_Insitu_Plasma",
        "anaesthetic": anaesthetic,
        "gwp_type": gwp_type,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)),
        "median": float(np.median(y)),
        "N": int(N),
    }
    
    if save_csv:
        filename = f"{outprefix}R11_sensitivity_{anaesthetic}_GWP{gwp_type}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")
    
    if make_plot:
        plot_file = f"{outprefix}R11_tornado_{anaesthetic}_GWP{gwp_type}.pdf"
        tornado_plot(df["Input"].tolist(), df["SRC"].values,
                    f"R11 Insitu Plasma Sensitivity ({anaesthetic.upper()}, GWP-{gwp_type})",
                    "Standardized Regression Coefficient (SRC)",
                    outfile=plot_file)
    
    return df, summary

def run_all_routes_sensitivity(seed=42, N=10000, anaesthetic='des', gwp_type='100',
                               save_csv=True, make_plot=True, outprefix="./"):
    """
    Run sensitivity analysis for all routes (R3-R11)
    """
    print("="*80)
    print(f"SENSITIVITY ANALYSIS: {anaesthetic.upper()}, GWP-{gwp_type}")
    print("="*80)
    
    results = {}
    
    print("\n--- Route R3: MSWI ---")
    df_r3, sum_r3 = run_R3_sensitivity(seed, N, anaesthetic, gwp_type, 
                                       save_csv, make_plot, outprefix)
    print(df_r3)
    print(f"\nOutput mean: {sum_r3['mean']:.2f} ± {sum_r3['std']:.2f} kgCO2e/kg")
    results['R3'] = (df_r3, sum_r3)
    
    print("\n--- Route R4: HCWI ---")
    df_r4, sum_r4 = run_R4_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r4)
    print(f"\nOutput mean: {sum_r4['mean']:.2f} ± {sum_r4['std']:.2f} kgCO2e/kg")
    results['R4'] = (df_r4, sum_r4)
    
    print("\n--- Route R5: Plasma ---")
    df_r5, sum_r5 = run_R5_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r5)
    print(f"\nOutput mean: {sum_r5['mean']:.2f} ± {sum_r5['std']:.2f} kgCO2e/kg")
    results['R5'] = (df_r5, sum_r5)
    
    print("\n--- Route R6: Continuous Release (Baseline) ---")
    df_r6, sum_r6 = run_R6_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    results['R6'] = (df_r6, sum_r6)
    
    print("\n--- Route R7: AC + MSWI ---")
    df_r7, sum_r7 = run_R7_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r7)
    print(f"\nOutput mean: {sum_r7['mean']:.2f} ± {sum_r7['std']:.2f} kgCO2e/kg")
    results['R7'] = (df_r7, sum_r7)
    
    print("\n--- Route R8: AC + HCWI ---")
    df_r8, sum_r8 = run_R8_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r8)
    print(f"\nOutput mean: {sum_r8['mean']:.2f} ± {sum_r8['std']:.2f} kgCO2e/kg")
    results['R8'] = (df_r8, sum_r8)
    
    print("\n--- Route R9: AC + Plasma ---")
    df_r9, sum_r9 = run_R9_sensitivity(seed, N, anaesthetic, gwp_type,
                                       save_csv, make_plot, outprefix)
    print(df_r9)
    print(f"\nOutput mean: {sum_r9['mean']:.2f} ± {sum_r9['std']:.2f} kgCO2e/kg")
    results['R9'] = (df_r9, sum_r9)
    
    print("\n--- Route R10: AC + On-site Incineration ---")
    df_r10, sum_r10 = run_R10_sensitivity(seed, N, anaesthetic, gwp_type,
                                         save_csv, make_plot, outprefix)
    print(df_r10)
    print(f"\nOutput mean: {sum_r10['mean']:.2f} ± {sum_r10['std']:.2f} kgCO2e/kg")
    results['R10'] = (df_r10, sum_r10)
    
    print("\n--- Route R11: In-situ Plasma ---")
    df_r11, sum_r11 = run_R11_sensitivity(seed, N, anaesthetic, gwp_type,
                                         save_csv, make_plot, outprefix)
    print(df_r11)
    print(f"\nOutput mean: {sum_r11['mean']:.2f} ± {sum_r11['std']:.2f} kgCO2e/kg")
    results['R11'] = (df_r11, sum_r11)
    
    print("\n" + "="*80)
    print("SUMMARY COMPARISON")
    print("="*80)
    print(f"{'Route':<25} {'Mean (kgCO2e/kg)':<20} {'Std Dev':<15}")
    print("-"*80)
    for route_name, (_, summary) in results.items():
        print(f"{route_name + ': ' + summary['route']:<25} {summary['mean']:<20.2f} {summary['std']:<15.2f}")
    print("="*80)
    
    return results

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
