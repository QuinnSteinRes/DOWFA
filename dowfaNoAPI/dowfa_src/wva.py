"""Waste volatile anaesthetic (WVA) pathway model.

Returns mean results for routes R6–R11.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np

from . import parameters as P
from .utils import rng, uniform_samples

@dataclass(frozen=True)
class WVAResult:
    total_years: int
    initial_population_m: float
    des_percent_mac: float
    n_samples: int
    totals_ktonnes_co2e: Dict[str, float]
    annual_ktonnes_co2e: Dict[str, float]
    perkg_kgco2e_per_kgVA: Dict[str, Dict[str, float]]  # route -> agent -> mean perkg
    stdev_perkg_kgco2e_per_kgVA: Dict[str, Dict[str, float]]  # route -> agent -> std perkg

def simulate_wva(
    total_years: int,
    initial_population_m: float,
    des_percent_mac: float,
    n_samples: int = 10_000,
    seed: Optional[int] = None,
) -> WVAResult:
    """Simulate WVA routes for a planning horizon.

    Parameters are taken from the original code. GWP fixed to 100y in the original.
    """
    if total_years < 0:
        raise ValueError("total_years must be non-negative")
    if initial_population_m < 0:
        raise ValueError("initial_population_m must be non-negative")
    if not (0 <= des_percent_mac <= 100):
        raise ValueError("des_percent_mac must be 0..100")

    r = rng(seed)
    cf4_gwp = P.cf4_100

    demunic = uniform_samples(r, P.min_demunic, P.max_demunic, n_samples)
    depharma = uniform_samples(r, P.min_depharma, P.max_depharma, n_samples)
    deplasma = uniform_samples(r, P.min_deplasma, P.max_deplasma, n_samples)

    munic_other = uniform_samples(r, P.min_munic_other, P.max_munic_other, n_samples)
    pharam_other = uniform_samples(r, P.min_pharam_other, P.max_pharam_other, n_samples)
    plasma_other = uniform_samples(r, P.min_plasma_other, P.max_plasma_other, n_samples)

    sur_100k = uniform_samples(r, P.min_sur_100k, P.max_sur_100k, n_samples)
    growth_rate = uniform_samples(r, P.min_growth_rate, P.max_growth_rate, n_samples)

    hgv_factor = uniform_samples(r, P.min_hgv_factor, P.max_hgv_factor, n_samples)
    ship_factor = uniform_samples(r, P.min_ship_factor, P.max_ship_factor, n_samples)

    network_distance = uniform_samples(r, P.min_network_distance, P.max_network_distance, n_samples)
    mwi_dist = uniform_samples(r, P.min_mwi_dist, P.max_mwi_dist, n_samples)
    pwi_dist = uniform_samples(r, P.min_pwi_dist, P.max_pwi_dist, n_samples)

    # 2017 baseline masses inferred from baseline ktCO2e and GWPs (as in original)
    kg_des = (P.ktco2e_des_base / P.des_100) * 1_000_000
    kg_sev = (P.ktco2e_sev_base / P.sev_100) * 1_000_000
    kg_iso = (P.ktco2e_iso_base / P.iso_100) * 1_000_000

    kmol_des = kg_des / P.des_gfm
    kmol_sev = kg_sev / P.sev_gfm
    kmol_iso = kg_iso / P.iso_gfm

    SI_des = kmol_des / P.MAC_des
    SI_sev = kmol_sev / P.MAC_sev
    SI_iso = kmol_iso / P.MAC_iso

    SI_baseline = (SI_des + SI_sev + SI_iso) / (P.sur_pro_base * 10)  # SI / 100k surgery

    ktco2e_per_SI_des = P.ktco2e_des_base / SI_des
    ktco2e_per_SI_sev = P.ktco2e_sev_base / SI_sev
    ktco2e_per_SI_iso = P.ktco2e_iso_base / SI_iso

    SI_100k_sur_iso = SI_iso / (P.sur_pro_base * 10)
    SI_100k_sur_des = (des_percent_mac/100.0) * SI_baseline
    SI_100k_sur_sev = SI_baseline - SI_100k_sur_des - SI_100k_sur_iso

    em_100k_des = ktco2e_per_SI_des * SI_100k_sur_des
    em_100k_sev = ktco2e_per_SI_sev * SI_100k_sur_sev
    em_100k_iso = ktco2e_per_SI_iso * SI_100k_sur_iso
    em_100k_total = em_100k_des + em_100k_sev + em_100k_iso

    # moles per 100k surgeries
    des_mol = SI_100k_sur_des * P.MAC_des
    sev_mol = SI_100k_sur_sev * P.MAC_sev
    iso_mol = SI_100k_sur_iso * P.MAC_iso

    total_mol = des_mol + sev_mol + iso_mol

    # Time loop (vectorized over samples)
    cumulative_surgeries = np.zeros(n_samples)
    for year in range(1, total_years + 1):
        pop_est = (initial_population_m * ((1.0 + growth_rate) ** year))
        surgeries_per_year = (pop_est / 100000.0) * sur_100k
        cumulative_surgeries += surgeries_per_year

    # masses of VA (kg) per 100k
    mass_va_des = des_mol * P.des_gfm
    mass_va_sev = sev_mol * P.sev_gfm
    mass_va_iso = iso_mol * P.iso_gfm
    mass_va_tonnes_100k = (mass_va_des + mass_va_sev + mass_va_iso) / 1000.0
    mass_va_tonnes = mass_va_tonnes_100k * (cumulative_surgeries * 10)

    # Activated carbon requirements
    tonnes_ac_req_100k = (des_mol / P.mol_per_kg_des) + (sev_mol / P.mol_per_kg_sev) + (iso_mol / P.mol_per_kg_iso)
    tonnes_ac_req = tonnes_ac_req_100k * (cumulative_surgeries * 10)

    perc_des_ac = ((des_mol / P.mol_per_kg_des) / tonnes_ac_req_100k) * 100.0
    perc_sev_ac = ((sev_mol / P.mol_per_kg_sev) / tonnes_ac_req_100k) * 100.0
    perc_iso_ac = ((iso_mol / P.mol_per_kg_iso) / tonnes_ac_req_100k) * 100.0

    tonnes_coconut = tonnes_ac_req * P.cococnut_to_ac_ratio

    ac_prod_em = P.ac_co2e_per_kg * (tonnes_ac_req / 1000.0)  # ktonnes CO2e

    ac_transport_1 = ((P.indo_to_sing_port * tonnes_coconut * ship_factor)/1e6) +                      ((P.sing_port_to_fac * tonnes_coconut * hgv_factor)/1e6) +                      ((P.sing_port_to_fac * tonnes_ac_req * hgv_factor)/1e6) +                      ((P.sing_to_eng_port * tonnes_ac_req * ship_factor)/1e6)

    ac_transport_2 = (network_distance * tonnes_ac_req * hgv_factor) / 1e6

    ac_transport_3_msw = (mwi_dist * (tonnes_ac_req + mass_va_tonnes) * hgv_factor) / 1e6
    ac_transport_3_pwi = (pwi_dist * (tonnes_ac_req + mass_va_tonnes) * hgv_factor) / 1e6
    ac_transport_3_pla = (network_distance * (tonnes_ac_req + mass_va_tonnes) * hgv_factor) / 1e6

    # Combustion proxy per 100k (vectorized over samples through de*)
    def combust_mol(mol, C, de):
        return (((mol * (1.0 - de)) * C * 44.1) / 1000.0) + ((((mol * de) * C * 88.0043) / 1000.0) * cf4_gwp)

    combust_des_msw = combust_mol(des_mol, P.des_atecon_C, demunic) * (cumulative_surgeries * 10)
    combust_des_pwi = combust_mol(des_mol, P.des_atecon_C, depharma) * (cumulative_surgeries * 10)
    combust_des_pla = combust_mol(des_mol, P.des_atecon_C, deplasma) * (cumulative_surgeries * 10)

    combust_sev_msw = combust_mol(sev_mol, P.sev_atecon_C, demunic) * (cumulative_surgeries * 10)
    combust_sev_pwi = combust_mol(sev_mol, P.sev_atecon_C, depharma) * (cumulative_surgeries * 10)
    combust_sev_pla = combust_mol(sev_mol, P.sev_atecon_C, deplasma) * (cumulative_surgeries * 10)

    combust_iso_msw = combust_mol(iso_mol, P.iso_atecon_C, demunic) * (cumulative_surgeries * 10)
    combust_iso_pwi = combust_mol(iso_mol, P.iso_atecon_C, depharma) * (cumulative_surgeries * 10)
    combust_iso_pla = combust_mol(iso_mol, P.iso_atecon_C, deplasma) * (cumulative_surgeries * 10)

    combust_msw_t = (combust_des_msw + combust_sev_msw + combust_iso_msw) / 1000.0
    combust_pwi_t = (combust_des_pwi + combust_sev_pwi + combust_iso_pwi) / 1000.0
    combust_pla_t = (combust_des_pla + combust_sev_pla + combust_iso_pla) / 1000.0

    em_per_kg_msw = (((mass_va_tonnes + tonnes_ac_req)*1000.0)*munic_other) / 1e6
    em_per_kg_pwi = (((mass_va_tonnes + tonnes_ac_req)*1000.0)*pharam_other) / 1e6
    em_per_kg_pla = (((mass_va_tonnes + tonnes_ac_req)*1000.0)*plasma_other) / 1e6

    # Route calculations (ktonnes CO2e)
    r6 = (cumulative_surgeries * 10) * em_100k_total
    r6_yr = r6 / total_years if total_years > 0 else r6*0.0

    r7 = ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_msw + combust_msw_t + em_per_kg_msw
    r7_yr = r7 / total_years if total_years > 0 else r7*0.0

    r8 = ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi + combust_pwi_t + em_per_kg_pwi
    r8_yr = r8 / total_years if total_years > 0 else r8*0.0

    r9 = ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pla + combust_pla_t + em_per_kg_pla
    r9_yr = r9 / total_years if total_years > 0 else r9*0.0

    r10 = ac_prod_em + ac_transport_1 + ac_transport_2 + combust_pwi_t + em_per_kg_pwi
    r10_yr = r10 / total_years if total_years > 0 else r10*0.0

    r11 = combust_pla_t + em_per_kg_pla
    r11_yr = r11 / total_years if total_years > 0 else r11*0.0

    totals = {
        "R6_release": float(np.mean(r6)),
        "R7_AC+MSWI": float(np.mean(r7)),
        "R8_AC+PWI": float(np.mean(r8)),
        "R9_AC+WM": float(np.mean(r9)),
        "R10_AC+onsite_inc": float(np.mean(r10)),
        "R11_in_situ_plasma": float(np.mean(r11)),
    }
    annual = {
        "R6_release": float(np.mean(r6_yr)) if total_years > 0 else 0.0,
        "R7_AC+MSWI": float(np.mean(r7_yr)) if total_years > 0 else 0.0,
        "R8_AC+PWI": float(np.mean(r8_yr)) if total_years > 0 else 0.0,
        "R9_AC+WM": float(np.mean(r9_yr)) if total_years > 0 else 0.0,
        "R10_AC+onsite_inc": float(np.mean(r10_yr)) if total_years > 0 else 0.0,
        "R11_in_situ_plasma": float(np.mean(r11_yr)) if total_years > 0 else 0.0,
    }

    # Per-kg per agent (kgCO2e/kgVA) — match original structure (route-level + base waste factor)
    # These are the same expressions as in your UI code but returned cleanly as means.
    # Note: mass_va_* are kg per 100k surgery.
    def perkg_route(route_total_ktonnes, combust_agent, perc_ac, other_factor, trans3, include_ac_prod_trans=True):
        # route_total_ktonnes not used; we keep your explicit per-agent algebra:
        pass

    # Original per-kg route values (directly from your code)
    r6_des = (em_100k_des / mass_va_des) * 1_000_000
    r6_sev = (em_100k_sev / mass_va_sev) * 1_000_000
    r6_iso = (em_100k_iso / mass_va_iso) * 1_000_000
    # Promote to arrays so std=0 is well-defined
    r6_des = np.full(n_samples, r6_des)
    r6_sev = np.full(n_samples, r6_sev)
    r6_iso = np.full(n_samples, r6_iso)

    r7_des = munic_other + (((((perc_des_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_des_msw)/1000.0)*1_000_000)/(mass_va_des*(cumulative_surgeries*10))
    r7_sev = munic_other + (((((perc_sev_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_sev_msw)/1000.0)*1_000_000)/(mass_va_sev*(cumulative_surgeries*10))
    r7_iso = munic_other + (((((perc_iso_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_iso_msw)/1000.0)*1_000_000)/(mass_va_iso*(cumulative_surgeries*10))

    r8_des = pharam_other + (((((perc_des_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_des_pwi)/1000.0)*1_000_000)/(mass_va_des*(cumulative_surgeries*10))
    r8_sev = pharam_other + (((((perc_sev_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_sev_pwi)/1000.0)*1_000_000)/(mass_va_sev*(cumulative_surgeries*10))
    r8_iso = pharam_other + (((((perc_iso_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_iso_pwi)/1000.0)*1_000_000)/(mass_va_iso*(cumulative_surgeries*10))

    r9_des = plasma_other + (((((perc_des_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_des_pla)/1000.0)*1_000_000)/(mass_va_des*(cumulative_surgeries*10))
    r9_sev = plasma_other + (((((perc_sev_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_sev_pla)/1000.0)*1_000_000)/(mass_va_sev*(cumulative_surgeries*10))
    r9_iso = plasma_other + (((((perc_iso_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_iso_pla)/1000.0)*1_000_000)/(mass_va_iso*(cumulative_surgeries*10))

    r10_des = pharam_other + (((((perc_des_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2))+combust_des_pwi)/1000.0)*1_000_000)/(mass_va_des*(cumulative_surgeries*10))
    r10_sev = pharam_other + (((((perc_sev_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2))+combust_sev_pwi)/1000.0)*1_000_000)/(mass_va_sev*(cumulative_surgeries*10))
    r10_iso = pharam_other + (((((perc_iso_ac/100.0)*(ac_prod_em + ac_transport_1 + ac_transport_2))+combust_iso_pwi)/1000.0)*1_000_000)/(mass_va_iso*(cumulative_surgeries*10))

    r11_des = plasma_other + (((combust_des_pla)/1000.0)*1_000_000)/(mass_va_des*(cumulative_surgeries*10))
    r11_sev = plasma_other + (((combust_sev_pla)/1000.0)*1_000_000)/(mass_va_sev*(cumulative_surgeries*10))
    r11_iso = plasma_other + (((combust_iso_pla)/1000.0)*1_000_000)/(mass_va_iso*(cumulative_surgeries*10))
    perkg = {
        "R6_release": {"des": float(np.mean(r6_des)), "sev": float(np.mean(r6_sev)), "iso": float(np.mean(r6_iso))},
        "R7_AC+MSWI": {"des": float(np.mean(r7_des)), "sev": float(np.mean(r7_sev)), "iso": float(np.mean(r7_iso))},
        "R8_AC+PWI": {"des": float(np.mean(r8_des)), "sev": float(np.mean(r8_sev)), "iso": float(np.mean(r8_iso))},
        "R9_AC+WM": {"des": float(np.mean(r9_des)), "sev": float(np.mean(r9_sev)), "iso": float(np.mean(r9_iso))},
        "R10_AC+onsite_inc": {"des": float(np.mean(r10_des)), "sev": float(np.mean(r10_sev)), "iso": float(np.mean(r10_iso))},
        "R11_in_situ_plasma": {"des": float(np.mean(r11_des)), "sev": float(np.mean(r11_sev)), "iso": float(np.mean(r11_iso))},
    }

    stdev_perkg = {
        "R6_release": {"des": float(np.std(r6_des)), "sev": float(np.std(r6_sev)), "iso": float(np.std(r6_iso))},
        "R7_AC+MSWI": {"des": float(np.std(r7_des)), "sev": float(np.std(r7_sev)), "iso": float(np.std(r7_iso))},
        "R8_AC+PWI": {"des": float(np.std(r8_des)), "sev": float(np.std(r8_sev)), "iso": float(np.std(r8_iso))},
        "R9_AC+WM": {"des": float(np.std(r9_des)), "sev": float(np.std(r9_sev)), "iso": float(np.std(r9_iso))},
        "R10_AC+onsite_inc": {"des": float(np.std(r10_des)), "sev": float(np.std(r10_sev)), "iso": float(np.std(r10_iso))},
        "R11_in_situ_plasma": {"des": float(np.std(r11_des)), "sev": float(np.std(r11_sev)), "iso": float(np.std(r11_iso))},
    }

    return WVAResult(
        total_years=int(total_years),
        initial_population_m=float(initial_population_m),
        des_percent_mac=float(des_percent_mac),
        n_samples=int(n_samples),
        totals_ktonnes_co2e=totals,
        annual_ktonnes_co2e=annual,
        perkg_kgco2e_per_kgVA=perkg,
        stdev_perkg_kgco2e_per_kgVA=stdev_perkg,
    )
