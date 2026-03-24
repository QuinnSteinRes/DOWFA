# distribution_comparison.py

import numpy as np
from parameters_m import *

# Number of samples for Monte Carlo
N_SAMPLES = 1000000

# ============================================================================
# DISTRIBUTION SAMPLING FUNCTIONS
# ============================================================================

def sample_parameter(min_val, max_val, distribution_type='uniform', n_samples=N_SAMPLES):
    """
    Sample a parameter according to specified distribution type.
    """
    
    if distribution_type == 'uniform':
        return np.random.uniform(min_val, max_val, n_samples)
    
    elif distribution_type == 'normal':
        mean = (max_val + min_val) / 2
        std = (max_val - min_val) / 6
        samples = np.random.normal(mean, std, n_samples)
        return np.clip(samples, min_val, max_val)
    
    elif distribution_type == 'lognormal':
        # For lognormal, parameters must be positive
        if min_val <= 0:
            raise ValueError(f"Lognormal distribution requires min_val > 0, got {min_val}")
        
        # Use method of moments to match desired mean and std to lognormal parameters
        desired_mean = (max_val + min_val) / 2
        desired_std = (max_val - min_val) / 6
        
        # Calculate lognormal parameters mu and sigma
        # E[X] = exp(mu + sigma^2/2) = desired_mean
        # Var[X] = (exp(sigma^2) - 1) * exp(2*mu + sigma^2) = desired_std^2
        
        variance = desired_std ** 2
        mu = np.log(desired_mean**2 / np.sqrt(variance + desired_mean**2))
        sigma = np.sqrt(np.log(variance / desired_mean**2 + 1))
        
        samples = np.random.lognormal(mu, sigma, n_samples)
        return np.clip(samples, min_val, max_val)
    
    elif distribution_type == 'triangular':
        mode = (max_val + min_val) / 2
        return np.random.triangular(min_val, mode, max_val, n_samples)
    
    else:
        raise ValueError(f"Unknown distribution type: {distribution_type}")


def sample_all_parameters(distribution_type='uniform', n_samples=N_SAMPLES):
    """
    Sample all uncertain parameters according to specified distribution.
    
    Returns dictionary with sampled arrays for each parameter.
    """
    
    samples = {}
    
    # Destruction efficiency parameters
    samples['demunic'] = sample_parameter(min_demunic, max_demunic, distribution_type, n_samples)
    samples['depharma'] = sample_parameter(min_depharma, max_depharma, distribution_type, n_samples)
    samples['deplasma'] = sample_parameter(min_deplasma, max_deplasma, distribution_type, n_samples)
    
    # General combustion emissions
    samples['munic_other'] = sample_parameter(min_munic_other, max_munic_other, distribution_type, n_samples)
    samples['pharam_other'] = sample_parameter(min_pharam_other, max_pharam_other, distribution_type, n_samples)
    samples['plasma_other'] = sample_parameter(min_plasma_other, max_plasma_other, distribution_type, n_samples)
    
    # Surgery and population parameters
    samples['sur_100k'] = sample_parameter(min_sur_100k, max_sur_100k, distribution_type, n_samples)
    samples['growth_rate'] = sample_parameter(min_growth_rate, max_growth_rate, distribution_type, n_samples)
    
    # Transportation factors
    samples['van_factor'] = sample_parameter(min_van_factor, max_van_factor, distribution_type, n_samples)
    samples['hgv_factor'] = sample_parameter(min_hgv_factor, max_hgv_factor, distribution_type, n_samples)
    samples['ship_factor'] = sample_parameter(min_ship_factor, max_ship_factor, distribution_type, n_samples)
    
    # Distances
    samples['network_distance'] = sample_parameter(min_network_distance, max_network_distance, distribution_type, n_samples)
    samples['mwi_dist'] = sample_parameter(min_mwi_dist, max_mwi_dist, distribution_type, n_samples)
    samples['pwi_dist'] = sample_parameter(min_pwi_dist, max_pwi_dist, distribution_type, n_samples)
    
    return samples


# ============================================================================
# LEGACY STOCK CALCULATIONS (R1-R5)
# ============================================================================

def calculate_legacy_stock_routes(samples, gwp_selection='GWP_100', variable_selection='des'):
    """
    Calculate emissions for legacy stock routes using provided parameter samples.
    
    Returns:
    --------
    dict : Dictionary with arrays for each route (R1-R5) in kgCO2e/kg
    """
    
    # GWP values
    gwp_categories = {
        "GWP_20": {"des": des_20, "sev": sev_20, "iso": iso_20},
        "GWP_100": {"des": des_100, "sev": sev_100, "iso": iso_100},
        "GWP_500": {"des": des_500, "sev": sev_500, "iso": iso_500}
    }
    
    cf4_gwp_values = {
        "GWP_20": cf4_20,
        "GWP_100": cf4_100,
        "GWP_500": cf4_500
    }
    
    # Volatile anaesthetic properties
    rho_values = {"des": rho_d, "sev": rho_s, "iso": rho_i}
    gfm_values = {"des": des_gfm, "sev": sev_gfm, "iso": iso_gfm}
    carbon_economy_values = {"des": des_atecon_C, "sev": sev_atecon_C, "iso": iso_atecon_C}
    
    # Get specific values
    rho = rho_values[variable_selection]
    gfm = gfm_values[variable_selection]
    carbon_economy = carbon_economy_values[variable_selection]
    cf4_gwp = cf4_gwp_values[gwp_selection]
    gwp_va = gwp_categories[gwp_selection][variable_selection]
    
    mol_val = 1000 / gfm
    
    # Combustion values (kgCO2e per kg VA)
    combust_value_msw = (
        (mol_val * (1 - samples['demunic']) * carbon_economy * 44.1) / 1000 +
        (mol_val * samples['demunic'] * carbon_economy * 88.0043 / 1000) * cf4_gwp
    )
    
    combust_value_pwi = (
        (mol_val * (1 - samples['depharma']) * carbon_economy * 44.1) / 1000 +
        (mol_val * samples['depharma'] * carbon_economy * 88.0043 / 1000) * cf4_gwp
    )
    
    combust_value_plasma = (
        (mol_val * (1 - samples['deplasma']) * carbon_economy * 44.1) / 1000 +
        (mol_val * samples['deplasma'] * carbon_economy * 88.0043 / 1000) * cf4_gwp
    )
    
    # Route calculations (per kg VA)
    results = {}
    
    # R1: Central storage (transport only)
    results['R1'] = (samples['network_distance'] * samples['van_factor']) / 1000
    
    # R2: Re-circulation (release)
    results['R2'] = np.full(len(samples['demunic']), gwp_va)
    
    # R3: Municipal waste incineration
    results['R3'] = (
        (samples['mwi_dist'] * samples['van_factor']) / 1000 +
        combust_value_msw +
        samples['munic_other']
    )
    
    # R4: Hazardous clinical waste incineration
    results['R4'] = (
        (samples['pwi_dist'] * samples['van_factor']) / 1000 +
        combust_value_pwi +
        samples['pharam_other']
    )
    
    # R5: Perfluorinated waste companies (plasma)
    results['R5'] = (
        (samples['network_distance'] * samples['van_factor']) / 1000 +
        combust_value_plasma +
        samples['plasma_other']
    )
    
    return results


# ============================================================================
# WASTE VOLATILE ANAESTHETIC CALCULATIONS (R6-R11)
# ============================================================================

def calculate_wva_routes(samples, initial_pop=66, total_years=26, des_prev=0.01, gwp_selection='GWP_100'):
    """
    Calculate emissions for waste volatile anaesthetic routes (R6-R11).
    
    Parameters:
    -----------
    samples : dict
        Dictionary of sampled parameter arrays
    initial_pop : float
        Initial population in millions
    total_years : int
        Number of years for projection
    des_prev : float
        Desflurane prevalence (as decimal, e.g., 0.01 for 1%)
    gwp_selection : str
        GWP time horizon ('GWP_100' typically)
        
    Returns:
    --------
    dict : Dictionary with results for each route and each anaesthetic
           Format: results[route][anaesthetic] = array of kgCO2e/kg values
    """
    
    cf4_gwp = cf4_100  # Using 100-year GWP
    
    # 2017 Baseline calculations
    kg_des = (ktco2e_des_base / des_100) * 1000000
    kg_sev = (ktco2e_sev_base / sev_100) * 1000000
    kg_iso = (ktco2e_iso_base / iso_100) * 1000000
    
    kmol_des = kg_des / des_gfm
    kmol_sev = kg_sev / sev_gfm
    kmol_iso = kg_iso / iso_gfm
    
    SI_des = kmol_des / MAC_des
    SI_sev = kmol_sev / MAC_sev
    SI_iso = kmol_iso / MAC_iso
    
    SI_baseline = (SI_des + SI_sev + SI_iso) / (sur_pro_base * 10)
    
    ktco2e_per_SI_des = ktco2e_des_base / SI_des
    ktco2e_per_SI_sev = ktco2e_sev_base / SI_sev
    ktco2e_per_SI_iso = ktco2e_iso_base / SI_iso
    
    # New emissions per 100k surgeries
    SI_100k_sur_iso = SI_iso / (sur_pro_base * 10)
    SI_100k_sur_des = (des_prev / 100) * SI_baseline
    SI_100k_sur_sev = SI_baseline - SI_100k_sur_des - SI_100k_sur_iso
    
    des_mol_calc = SI_100k_sur_des * MAC_des
    sev_mol_calc = SI_100k_sur_sev * MAC_sev
    iso_mol_calc = SI_100k_sur_iso * MAC_iso
    
    # Calculate cumulative surgeries over time period
    cumulative_surgeries = np.zeros_like(samples['growth_rate'])
    for year in range(1, total_years + 1):
        population_estimate = initial_pop * ((1 + samples['growth_rate']) ** year)
        surgeries_per_year = (population_estimate / 100000) * samples['sur_100k']
        cumulative_surgeries += surgeries_per_year
    
    # Activated Carbon costs (per 100k surgeries)
    mass_va_des = des_mol_calc * des_gfm
    mass_va_sev = sev_mol_calc * sev_gfm
    mass_va_iso = iso_mol_calc * iso_gfm
    mass_va_tonnes_100k = (mass_va_des + mass_va_sev + mass_va_iso) / 1000
    mass_va_tonnes = mass_va_tonnes_100k * (cumulative_surgeries * 10)
    
    tonnes_ac_req_100k = ((des_mol_calc / mol_per_kg_des) + 
                          (sev_mol_calc / mol_per_kg_sev) + 
                          (iso_mol_calc / mol_per_kg_iso))
    tonnes_ac_req = tonnes_ac_req_100k * (cumulative_surgeries * 10)
    
    perc_des_ac = ((des_mol_calc / mol_per_kg_des) / tonnes_ac_req_100k) * 100
    perc_sev_ac = ((sev_mol_calc / mol_per_kg_sev) / tonnes_ac_req_100k) * 100
    perc_iso_ac = ((iso_mol_calc / mol_per_kg_iso) / tonnes_ac_req_100k) * 100
    
    tonnnes_coconut = tonnes_ac_req * cococnut_to_ac_ratio
    
    # AC supply chain emissions (ktonnesCO2e)
    ac_production_emmisions = ac_co2e_per_kg * (tonnes_ac_req / 1000)
    
    ac_transport_1 = (
        ((indo_to_sing_port * tonnnes_coconut * samples['ship_factor']) / 1000000) +
        ((sing_port_to_fac * tonnnes_coconut * samples['hgv_factor']) / 1000000) +
        ((sing_port_to_fac * tonnes_ac_req * samples['hgv_factor']) / 1000000) +
        ((sing_to_eng_port * tonnes_ac_req * samples['ship_factor']) / 1000000)
    )
    
    ac_transport_2 = (samples['network_distance'] * tonnes_ac_req * samples['hgv_factor']) / 1000000
    
    ac_transport_3_msw = (samples['mwi_dist'] * (tonnes_ac_req + mass_va_tonnes) * samples['hgv_factor']) / 1000000
    ac_transport_3_pwi = (samples['pwi_dist'] * (tonnes_ac_req + mass_va_tonnes) * samples['hgv_factor']) / 1000000
    ac_transport_3_pla = (samples['network_distance'] * (tonnes_ac_req + mass_va_tonnes) * samples['hgv_factor']) / 1000000
    
    # Combustion emissions for each anaesthetic type
    # Desflurane
    combust_value_des_msw = (
        (((des_mol_calc * (1 - samples['demunic']) * des_atecon_C * 44.1)) / 1000) +
        (((((des_mol_calc * (samples['demunic'])) * des_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_des_msw = combust_value_des_msw * (cumulative_surgeries * 10)
    
    combust_value_des_pwi = (
        (((des_mol_calc * (1 - samples['depharma']) * des_atecon_C * 44.1)) / 1000) +
        (((((des_mol_calc * (samples['depharma'])) * des_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_des_pwi = combust_value_des_pwi * (cumulative_surgeries * 10)
    
    combust_value_des_pla = (
        (((des_mol_calc * (1 - samples['deplasma']) * des_atecon_C * 44.1)) / 1000) +
        (((((des_mol_calc * (samples['deplasma'])) * des_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_des_pla = combust_value_des_pla * (cumulative_surgeries * 10)
    
    # Sevoflurane
    combust_value_sev_msw = (
        (((sev_mol_calc * (1 - samples['demunic']) * sev_atecon_C * 44.1)) / 1000) +
        (((((sev_mol_calc * (samples['demunic'])) * sev_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_sev_msw = combust_value_sev_msw * (cumulative_surgeries * 10)
    
    combust_value_sev_pwi = (
        (((sev_mol_calc * (1 - samples['depharma']) * sev_atecon_C * 44.1)) / 1000) +
        (((((sev_mol_calc * (samples['depharma'])) * sev_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_sev_pwi = combust_value_sev_pwi * (cumulative_surgeries * 10)
    
    combust_value_sev_pla = (
        (((sev_mol_calc * (1 - samples['deplasma']) * sev_atecon_C * 44.1)) / 1000) +
        (((((sev_mol_calc * (samples['deplasma'])) * sev_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_sev_pla = combust_value_sev_pla * (cumulative_surgeries * 10)
    
    # Isoflurane
    combust_value_iso_msw = (
        (((iso_mol_calc * (1 - samples['demunic']) * iso_atecon_C * 44.1)) / 1000) +
        (((((iso_mol_calc * (samples['demunic'])) * iso_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_iso_msw = combust_value_iso_msw * (cumulative_surgeries * 10)
    
    combust_value_iso_pwi = (
        (((iso_mol_calc * (1 - samples['depharma']) * iso_atecon_C * 44.1)) / 1000) +
        (((((iso_mol_calc * (samples['depharma'])) * iso_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_iso_pwi = combust_value_iso_pwi * (cumulative_surgeries * 10)
    
    combust_value_iso_pla = (
        (((iso_mol_calc * (1 - samples['deplasma']) * iso_atecon_C * 44.1)) / 1000) +
        (((((iso_mol_calc * (samples['deplasma'])) * iso_atecon_C * 88.0043)) / 1000) * cf4_gwp)
    )
    combust_emmisions_iso_pla = combust_value_iso_pla * (cumulative_surgeries * 10)
    
    # Emissions from kg waste incineration
    emmisions_per_kg_msw = (((mass_va_tonnes + tonnes_ac_req) * 1000) * samples['munic_other']) / 1000000
    emmisions_per_kg_pwi = (((mass_va_tonnes + tonnes_ac_req) * 1000) * samples['pharam_other']) / 1000000
    emmisions_per_kg_pla = (((mass_va_tonnes + tonnes_ac_req) * 1000) * samples['plasma_other']) / 1000000
    
    # Calculate per kg VA emissions for each route
    results = {}
    
    # R6: Continuous Release (baseline for each anaesthetic)
    results['R6'] = {
        'des': np.full_like(samples['demunic'], des_100),
        'sev': np.full_like(samples['demunic'], sev_100),
        'iso': np.full_like(samples['demunic'], iso_100)
    }
    
    # R7: AC capture with municipal waste incineration
    results['R7'] = {
        'des': samples['munic_other'] + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_emmisions_des_msw)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10))),
        'sev': samples['munic_other'] + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_emmisions_sev_msw)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10))),
        'iso': samples['munic_other'] + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_emmisions_iso_msw)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
    }
    
    # R8: AC capture with hazardous clinical waste incineration
    results['R8'] = {
        'des': samples['pharam_other'] + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_emmisions_des_pwi)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10))),
        'sev': samples['pharam_other'] + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_emmisions_sev_pwi)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10))),
        'iso': samples['pharam_other'] + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_emmisions_iso_pwi)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
    }
    
    # R9: AC capture with perfluorinated waste companies
    results['R9'] = {
        'des': samples['plasma_other'] + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_emmisions_des_pla)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10))),
        'sev': samples['plasma_other'] + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_emmisions_sev_pla)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10))),
        'iso': samples['plasma_other'] + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_emmisions_iso_pla)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
    }
    
    # R10: AC capture with on-site incineration
    results['R10'] = {
        'des': samples['pharam_other'] + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2))+combust_emmisions_des_pwi)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10))),
        'sev': samples['pharam_other'] + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2))+combust_emmisions_sev_pwi)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10))),
        'iso': samples['pharam_other'] + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2))+combust_emmisions_iso_pwi)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
    }
    
    # R11: In-situ plasma destruction
    results['R11'] = {
        'des': samples['plasma_other'] + ((((combust_emmisions_des_pla)/1000)*1000000)/((mass_va_des * (cumulative_surgeries*10)))),
        'sev': samples['plasma_other'] + ((((combust_emmisions_sev_pla)/1000)*1000000)/((mass_va_sev * (cumulative_surgeries*10)))),
        'iso': samples['plasma_other'] + ((((combust_emmisions_iso_pla)/1000)*1000000)/((mass_va_iso * (cumulative_surgeries*10))))
    }
    
    return results


# ============================================================================
# COMPARISON ANALYSIS
# ============================================================================

def compare_legacy_stock_distributions(routes_dict, distribution_names):
    """
    Compare legacy stock route results across different distributions.
    """
    
    # Get route names from first distribution
    first_dist = distribution_names[0]
    route_names = list(routes_dict[first_dist].keys())
    
    print("\n" + "="*100)
    print("DISTRIBUTION COMPARISON: LEGACY STOCK ROUTES (R1-R5) - Desflurane, GWP_100")
    print("="*100)
    print(f"{'Route':<10}", end='')
    
    for dist_name in distribution_names:
        print(f"{dist_name.upper():<35}", end='')
    print()
    
    print(f"{'':10}", end='')
    for _ in distribution_names:
        print(f"{'Mean ± Std (kgCO2e/kg)':<35}", end='')
    print()
    
    print("-"*100)
    
    # Calculate and display statistics for each route
    for route in route_names:
        print(f"{route:<10}", end='')
        
        route_stats = {}
        for dist_name in distribution_names:
            values = routes_dict[dist_name][route]
            mean = np.mean(values)
            std = np.std(values)
            route_stats[dist_name] = (mean, std)
            print(f"{mean:>12.2f} ± {std:<18.2f}", end='')
        
        print()
    
    print("="*100)
    
    # Show percentage differences from uniform (baseline)
    if len(distribution_names) > 1:
        print("\n" + "="*100)
        print("PERCENTAGE DIFFERENCE FROM UNIFORM DISTRIBUTION")
        print("="*100)
        print(f"{'Route':<10}", end='')
        
        for dist_name in distribution_names[1:]:  # Skip uniform
            print(f"{dist_name.upper():<35}", end='')
        print()
        
        print(f"{'':10}", end='')
        for _ in distribution_names[1:]:
            print(f"{'% Diff (Mean) | % Diff (Std)':<35}", end='')
        print()
        
        print("-"*100)
        
        for route in route_names:
            print(f"{route:<10}", end='')
            
            uniform_mean = np.mean(routes_dict['uniform'][route])
            uniform_std = np.std(routes_dict['uniform'][route])
            
            for dist_name in distribution_names[1:]:
                test_mean = np.mean(routes_dict[dist_name][route])
                test_std = np.std(routes_dict[dist_name][route])
                
                # Calculate percentage differences
                if uniform_mean != 0:
                    mean_diff = ((test_mean - uniform_mean) / uniform_mean) * 100
                else:
                    mean_diff = 0.0
                    
                if uniform_std != 0:
                    std_diff = ((test_std - uniform_std) / uniform_std) * 100
                else:
                    std_diff = 0.0
                
                print(f"{mean_diff:>10.2f}%  | {std_diff:>10.2f}%       ", end='')
            
            print()
        
        print("="*100)


def compare_wva_distributions(routes_dict, distribution_names):
    """
    Compare WVA route results across different distributions.
    """
    
    # Get route names from first distribution
    first_dist = distribution_names[0]
    route_names = list(routes_dict[first_dist].keys())
    anaesthetics = ['des', 'sev', 'iso']
    
    for anaesthetic in anaesthetics:
        print("\n" + "="*120)
        print(f"DISTRIBUTION COMPARISON: WVA ROUTES (R6-R11) - {anaesthetic.upper()}, GWP_100")
        print("="*120)
        print(f"{'Route':<10}", end='')
        
        for dist_name in distribution_names:
            print(f"{dist_name.upper():<40}", end='')
        print()
        
        print(f"{'':10}", end='')
        for _ in distribution_names:
            print(f"{'Mean ± Std (kgCO2e/kg)':<40}", end='')
        print()
        
        print("-"*120)
        
        # Calculate and display statistics for each route
        for route in route_names:
            print(f"{route:<10}", end='')
            
            for dist_name in distribution_names:
                values = routes_dict[dist_name][route][anaesthetic]
                mean = np.mean(values)
                std = np.std(values)
                print(f"{mean:>15.2f} ± {std:<22.2f}", end='')
            
            print()
        
        print("="*120)
        
        # Show percentage differences from uniform (baseline)
        if len(distribution_names) > 1:
            print("\n" + "="*120)
            print("PERCENTAGE DIFFERENCE FROM UNIFORM DISTRIBUTION")
            print("="*120)
            print(f"{'Route':<10}", end='')
            
            for dist_name in distribution_names[1:]:  # Skip uniform
                print(f"{dist_name.upper():<40}", end='')
            print()
            
            print(f"{'':10}", end='')
            for _ in distribution_names[1:]:
                print(f"{'% Diff (Mean) | % Diff (Std)':<40}", end='')
            print()
            
            print("-"*120)
            
            for route in route_names:
                print(f"{route:<10}", end='')
                
                uniform_mean = np.mean(routes_dict['uniform'][route][anaesthetic])
                uniform_std = np.std(routes_dict['uniform'][route][anaesthetic])
                
                for dist_name in distribution_names[1:]:
                    test_mean = np.mean(routes_dict[dist_name][route][anaesthetic])
                    test_std = np.std(routes_dict[dist_name][route][anaesthetic])
                    
                    # Calculate percentage differences
                    if uniform_mean != 0:
                        mean_diff = ((test_mean - uniform_mean) / uniform_mean) * 100
                    else:
                        mean_diff = 0.0
                        
                    if uniform_std != 0:
                        std_diff = ((test_std - uniform_std) / uniform_std) * 100
                    else:
                        std_diff = 0.0
                    
                    print(f"{mean_diff:>12.2f}%  | {std_diff:>12.2f}%         ", end='')
                
                print()
            
            print("="*120)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    print("\n" + "="*100)
    print("DISTRIBUTION COMPARISON ANALYSIS - ALL ROUTES")
    print("="*100)
    print(f"Number of Monte Carlo samples: {N_SAMPLES}")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Distribution types to compare
    distribution_types = ['uniform', 'normal', 'lognormal']
    
    # ========================================================================
    # LEGACY STOCK ROUTES (R1-R5)
    # ========================================================================
    
    print("\n" + "-"*100)
    print("PART 1: LEGACY STOCK ROUTES (R1-R5)")
    print("-"*100)
    
    legacy_results = {}
    
    for dist_type in distribution_types:
        print(f"\nSampling parameters using {dist_type} distribution...")
        samples = sample_all_parameters(distribution_type=dist_type, n_samples=N_SAMPLES)
        
        print(f"Calculating legacy stock routes...")
        routes = calculate_legacy_stock_routes(samples, gwp_selection='GWP_100', variable_selection='des')
        
        legacy_results[dist_type] = routes
    
    # Display comparison
    compare_legacy_stock_distributions(legacy_results, distribution_types)
    
    # ========================================================================
    # WASTE VOLATILE ANAESTHETIC ROUTES (R6-R11)
    # ========================================================================
    
    print("\n" + "-"*100)
    print("PART 2: WASTE VOLATILE ANAESTHETIC ROUTES (R6-R11)")
    print("-"*100)
    
    wva_results = {}
    
    for dist_type in distribution_types:
        print(f"\nSampling parameters using {dist_type} distribution...")
        samples = sample_all_parameters(distribution_type=dist_type, n_samples=N_SAMPLES)
        
        print(f"Calculating WVA routes...")
        routes = calculate_wva_routes(samples, initial_pop=66, total_years=26, des_prev=0.01)
        
        wva_results[dist_type] = routes
    
    # Display comparison
    compare_wva_distributions(wva_results, distribution_types)
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    
    
    print("\n✓ Analysis complete for all routes (R1-R11)!")