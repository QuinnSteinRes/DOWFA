"""Build emissions time series for WVA disposal routes (R6–R11) for temperature modelling.

This module mirrors the parameterisation already present in `dowfa.wva` (the CO2e model),
but returns species-resolved emissions time series (CO2 / DES / SEV / ISO / CF4) suitable for
the temperature impulse-response model.

Key idea
--------
`dowfa.wva` infers a *baseline* anaesthetic activity level from published ktCO2e baselines
and agent GWPs, converts this to a "surgery index" per 100k surgeries, then projects activity
forward using:
- `sur_100k`  : surgeries per 100k population per year (uncertain; sampled)
- `growth_rate`: annual population growth (uncertain; sampled)

We use the same structure, but compute *annual* agent masses (kg/year) and then translate each
disposal route into time series of species emissions.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Tuple
import numpy as np

from . import parameters as P
from .utils import rng, uniform_samples

@dataclass(frozen=True)
class WVAEmissions:
    t_years: np.ndarray
    emissions_g_per_year: Dict[str, np.ndarray]  # species -> rate (g/yr)

def _annual_va_mass_kg_per_agent(
    total_years: int,
    initial_population_m: float,
    des_percent_mac: float,
    sur_100k: np.ndarray,
    growth_rate: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute annual VA mass trajectories (kg/year) for DES/SEV/ISO.

    This follows the same baseline derivation as `dowfa.wva`, but outputs per-year mass.
    Shapes:
      - sur_100k, growth_rate: (n_samples,)
      - returns: kg_des, kg_sev, kg_iso of shape (n_samples, total_years+1) with year 0 included.
    """
    n_samples = len(sur_100k)
    years = np.arange(0, total_years + 1, dtype=float)  # include year 0

    # Baseline masses inferred from baseline ktCO2e and GWPs (as in dowfa.wva)
    kg_des_base = (P.ktco2e_des_base / P.des_100) * 1_000_000
    kg_sev_base = (P.ktco2e_sev_base / P.sev_100) * 1_000_000
    kg_iso_base = (P.ktco2e_iso_base / P.iso_100) * 1_000_000

    kmol_des = kg_des_base / P.des_gfm
    kmol_sev = kg_sev_base / P.sev_gfm
    kmol_iso = kg_iso_base / P.iso_gfm

    SI_des = kmol_des / P.MAC_des
    SI_sev = kmol_sev / P.MAC_sev
    SI_iso = kmol_iso / P.MAC_iso

    SI_baseline = (SI_des + SI_sev + SI_iso) / (P.sur_pro_base * 10)  # SI per 100k surgeries

    SI_100k_sur_iso = SI_iso / (P.sur_pro_base * 10)
    SI_100k_sur_des = (des_percent_mac / 100.0) * SI_baseline
    SI_100k_sur_sev = SI_baseline - SI_100k_sur_des - SI_100k_sur_iso

    # moles per 100k surgeries (kmol if MAC_* is kmol/SI)
    des_mol_100k = SI_100k_sur_des * P.MAC_des
    sev_mol_100k = SI_100k_sur_sev * P.MAC_sev
    iso_mol_100k = SI_100k_sur_iso * P.MAC_iso

    # mass per 100k surgeries (kg)
    mass_des_100k = des_mol_100k * P.des_gfm
    mass_sev_100k = sev_mol_100k * P.sev_gfm
    mass_iso_100k = iso_mol_100k * P.iso_gfm

    # Annual surgeries (count) for each sample and year: surgeries/year = (pop/100k)*sur_100k
    pop0 = initial_population_m  # in millions
    pop_est_m = pop0 * (1.0 + growth_rate[:, None]) ** years[None, :]
    surgeries_per_year = (pop_est_m * 1e6 / 100000.0) * sur_100k[:, None]  # count/year

    # Convert surgeries count to "per 100k surgeries" multiplier
    mult_100k = surgeries_per_year / 100000.0

    kg_des = mass_des_100k * mult_100k
    kg_sev = mass_sev_100k * mult_100k
    kg_iso = mass_iso_100k * mult_100k

    return kg_des, kg_sev, kg_iso

def wva_route_emissions(
    route: str,
    total_years: int = 50,
    initial_population_m: float = 45.0,
    des_percent_mac: float = 12.0,
    n_samples: int = 1000,
    seed: Optional[int] = None,
    dt_years: float = 1.0,
    cf4_yield_fraction: Optional[float] = None,
) -> Dict[str, WVAEmissions]:
    """Generate Monte Carlo emissions schedules for a given WVA route (R6–R11)."""
    route = route.strip()
    if total_years <= 0:
        raise ValueError("total_years must be positive")

    r = rng(seed)

    # Uncertain drivers (same ranges as wva.py)
    sur_100k = uniform_samples(r, P.min_sur_100k, P.max_sur_100k, n_samples)
    growth_rate = uniform_samples(r, P.min_growth_rate, P.max_growth_rate, n_samples)

    # Other uncertainties (transport + destruction proxies)
    van_factor = uniform_samples(r, P.min_van_factor, P.max_van_factor, n_samples)
    mwi_dist = uniform_samples(r, P.min_mwi_dist, P.max_mwi_dist, n_samples)
    pwi_dist = uniform_samples(r, P.min_pwi_dist, P.max_pwi_dist, n_samples)
    network_distance = uniform_samples(r, P.min_network_distance, P.max_network_distance, n_samples)

    demunic = uniform_samples(r, P.min_demunic, P.max_demunic, n_samples)
    depharma = uniform_samples(r, P.min_depharma, P.max_depharma, n_samples)
    deplasma = uniform_samples(r, P.min_deplasma, P.max_deplasma, n_samples)

    munic_other = uniform_samples(r, P.min_munic_other, P.max_munic_other, n_samples)
    pharam_other = uniform_samples(r, P.min_pharam_other, P.max_pharam_other, n_samples)
    plasma_other = uniform_samples(r, P.min_plasma_other, P.max_plasma_other, n_samples)

    # Annual agent masses (kg/year), include year 0..total_years
    kg_des, kg_sev, kg_iso = _annual_va_mass_kg_per_agent(
        total_years=total_years,
        initial_population_m=initial_population_m,
        des_percent_mac=des_percent_mac,
        sur_100k=sur_100k,
        growth_rate=growth_rate,
    )
    t = np.arange(0, total_years + 1, dtype=float) * dt_years
    n_t = len(t)

    # Conversion helpers
    gfm = {"des": P.des_gfm, "sev": P.sev_gfm, "iso": P.iso_gfm}
    atecon_C = {"des": P.des_atecon_C, "sev": P.sev_atecon_C, "iso": P.iso_atecon_C}
    MW_CO2 = 44.01
    MW_CF4 = 88.00

    def combustion_products(agent: str, de: np.ndarray, kg_series: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        de_eff = de if cf4_yield_fraction is None else np.full_like(de, cf4_yield_fraction)
        mol = (kg_series * 1000.0) / gfm[agent]  # mol/year, shape (n_samples, t)
        g_co2 = mol * (1.0 - de_eff[:, None]) * atecon_C[agent] * MW_CO2
        g_cf4 = mol * de_eff[:, None] * atecon_C[agent] * MW_CF4
        return g_co2, g_cf4

    def co2_transport(dist_km: np.ndarray, kg_total_series: np.ndarray) -> np.ndarray:
        tonneVA = (kg_total_series / 1000.0)
        kgco2 = dist_km[:, None] * tonneVA * van_factor[:, None]
        return kgco2 * 1000.0

    def co2_other(kg_total_series: np.ndarray, other_factor: np.ndarray) -> np.ndarray:
        kgco2 = kg_total_series * other_factor[:, None]
        return kgco2 * 1000.0

    out: Dict[str, WVAEmissions] = {f"sample_{i:05d}": WVAEmissions(t_years=t, emissions_g_per_year={}) for i in range(n_samples)}

    def set_species(spec: str, arr: np.ndarray):
        for i in range(n_samples):
            out[f"sample_{i:05d}"].emissions_g_per_year[spec] = arr[i]

    # Base agent emissions (g/year)
    DES_g = kg_des * 1000.0
    SEV_g = kg_sev * 1000.0
    ISO_g = kg_iso * 1000.0
    KG_total = kg_des + kg_sev + kg_iso

    if route == "R6_release":
        set_species("DES", DES_g)
        set_species("SEV", SEV_g)
        set_species("ISO", ISO_g)
        return out

    if route == "R7_AC+MSWI":
        gco2 = co2_transport(mwi_dist, KG_total) + co2_other(KG_total, munic_other)
        gco2_des, gcf4_des = combustion_products("des", demunic, kg_des)
        gco2_sev, gcf4_sev = combustion_products("sev", demunic, kg_sev)
        gco2_iso, gcf4_iso = combustion_products("iso", demunic, kg_iso)
        set_species("CO2", gco2 + gco2_des + gco2_sev + gco2_iso)
        set_species("CF4", gcf4_des + gcf4_sev + gcf4_iso)
        return out

    if route == "R8_AC+PWI":
        gco2 = co2_transport(pwi_dist, KG_total) + co2_other(KG_total, pharam_other)
        gco2_des, gcf4_des = combustion_products("des", depharma, kg_des)
        gco2_sev, gcf4_sev = combustion_products("sev", depharma, kg_sev)
        gco2_iso, gcf4_iso = combustion_products("iso", depharma, kg_iso)
        set_species("CO2", gco2 + gco2_des + gco2_sev + gco2_iso)
        set_species("CF4", gcf4_des + gcf4_sev + gcf4_iso)
        return out

    if route == "R9_AC+WM":
        gco2 = co2_transport(network_distance, KG_total) + co2_other(KG_total, plasma_other)
        gco2_des, gcf4_des = combustion_products("des", deplasma, kg_des)
        gco2_sev, gcf4_sev = combustion_products("sev", deplasma, kg_sev)
        gco2_iso, gcf4_iso = combustion_products("iso", deplasma, kg_iso)
        set_species("CO2", gco2 + gco2_des + gco2_sev + gco2_iso)
        set_species("CF4", gcf4_des + gcf4_sev + gcf4_iso)
        return out

    if route == "R10_AC+onsite_inc":
        return wva_route_emissions("R8_AC+PWI", total_years, initial_population_m, des_percent_mac, n_samples, seed, dt_years, cf4_yield_fraction)

    if route == "R11_in_situ_plasma":
        return wva_route_emissions("R9_AC+WM", total_years, initial_population_m, des_percent_mac, n_samples, seed, dt_years, cf4_yield_fraction)

    raise ValueError(f"Unknown route: {route}")
