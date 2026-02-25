"""Build emissions time series for legacy stock routes for use with the temperature model."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np

from . import parameters as P
from .utils import rng, uniform_samples, agent_key, gwp_horizon_key

@dataclass(frozen=True)
class LegacyEmissions:
    t_years: np.ndarray
    emissions_g_per_year: Dict[str, np.ndarray]

def legacy_route_emissions_pulse(
    litres: float,
    agent: str = "des",
    route: str = "R2_recirculation",
    n_samples: int = 1000,
    seed: Optional[int] = None,
    years_total: float = 200.0,
    dt_years: float = 0.25,
    cf4_yield_fraction: Optional[float] = None,
) -> Dict[str, LegacyEmissions]:
    """Generate per-route emissions schedules for legacy stock disposal, with Monte Carlo.

    For temperature modelling we need emissions by *species* (CO2, VA, CF4 ...). We adopt:
    - Transport and 'other emissions' terms are treated as CO2 mass (g).
    - Recirculation/release: emitted as parent VA mass.
    - Incineration/plasma: parent VA is destroyed; combustion proxy produces CO2 and CF4.
      The original model uses a proxy 'de' that scales CF4 formation; we keep that here.
    - Optionally, override CF4 formation yield by providing cf4_yield_fraction (mass fraction
      of carbon converted to CF4 proxy), for sensitivity studies.

    Returns a dict mapping sample label -> LegacyEmissions.
    """
    if litres < 0:
        raise ValueError("litres must be non-negative")

    a = agent_key(agent)
    route = route.strip()

    rho_values = {"des": P.rho_d, "sev": P.rho_s, "iso": P.rho_i}
    gfm_values = {"des": P.des_gfm, "sev": P.sev_gfm, "iso": P.iso_gfm}
    carbon_economy_values = {"des": P.des_atecon_C, "sev": P.sev_atecon_C, "iso": P.iso_atecon_C}

    rho = rho_values[a]
    gfm = gfm_values[a]
    carbon_economy = carbon_economy_values[a]

    kg_va = litres * rho / 1000.0
    g_va = kg_va * 1000.0

    mol_per_kg = 1000.0 / gfm  # mol/kg

    r = rng(seed)
    van_factor = uniform_samples(r, P.min_van_factor, P.max_van_factor, n_samples)
    network_distance = uniform_samples(r, P.min_network_distance, P.max_network_distance, n_samples)
    mwi_dist = uniform_samples(r, P.min_mwi_dist, P.max_mwi_dist, n_samples)
    pwi_dist = uniform_samples(r, P.min_pwi_dist, P.max_pwi_dist, n_samples)

    demunic = uniform_samples(r, P.min_demunic, P.max_demunic, n_samples)
    depharma = uniform_samples(r, P.min_depharma, P.max_depharma, n_samples)
    deplasma = uniform_samples(r, P.min_deplasma, P.max_deplasma, n_samples)

    munic_other = uniform_samples(r, P.min_munic_other, P.max_munic_other, n_samples)
    pharam_other = uniform_samples(r, P.min_pharam_other, P.max_pharam_other, n_samples)
    plasma_other = uniform_samples(r, P.min_plasma_other, P.max_plasma_other, n_samples)

    # Time grid
    t = np.arange(int(np.round(years_total/dt_years)) + 1) * dt_years
    n_t = len(t)

    def pulse_series(total_g: np.ndarray) -> np.ndarray:
        # Convert pulse mass (g) at t=0 into rate over first time step (g/yr)
        rate = np.zeros((n_samples, n_t))
        rate[:, 0] = total_g / dt_years
        return rate

    # R1 central storage: CO2 only (transport proxy): x1*x2/1000 kgCO2/kg -> convert to gCO2
    r1_perkg = (network_distance * van_factor) / 1000.0  # kgCO2/kgVA
    r1_g_co2 = (r1_perkg * kg_va) * 1000.0

    # R2 recirculation: emit parent VA
    r2_g_va = np.full(n_samples, g_va)

    # Combustion proxy: produce CO2 and CF4 using de
    # CO2 mass (g) from carbon economy * (1-de) * mol * MW_CO2
    # CF4 mass (g) from carbon economy * (de) * mol * MW_CF4
    MW_CO2 = 44.01
    MW_CF4 = 88.00
    mol_total = mol_per_kg * kg_va  # mol

    def combustion_products(de: np.ndarray):
        de_eff = de if cf4_yield_fraction is None else np.full_like(de, cf4_yield_fraction)
        g_co2 = (mol_total * (1.0 - de_eff) * carbon_economy * MW_CO2)
        g_cf4 = (mol_total * de_eff * carbon_economy * MW_CF4)
        return g_co2, g_cf4

    # R3 municipal incineration: transport CO2 + combustion CO2/CF4 + other as CO2
    gco2_comb3, gcf4_3 = combustion_products(demunic)
    gco2_trans3 = (mwi_dist * (kg_va/1000.0) * van_factor) * 1000.0  # kg->g
    gco2_other3 = (kg_va * munic_other) * 1000.0
    r3_g_co2 = gco2_trans3 + gco2_other3 + gco2_comb3
    r3_g_cf4 = gcf4_3

    # R4 hazardous incineration: transport + combustion + other
    gco2_comb4, gcf4_4 = combustion_products(depharma)
    gco2_trans4 = (pwi_dist * (kg_va/1000.0) * van_factor) * 1000.0
    gco2_other4 = (kg_va * pharam_other) * 1000.0
    r4_g_co2 = gco2_trans4 + gco2_other4 + gco2_comb4
    r4_g_cf4 = gcf4_4

    # R5 perfluoro company/plasma: network transport + combustion + other
    gco2_comb5, gcf4_5 = combustion_products(deplasma)
    gco2_trans5 = (network_distance * (kg_va/1000.0) * van_factor) * 1000.0
    gco2_other5 = (kg_va * plasma_other) * 1000.0
    r5_g_co2 = gco2_trans5 + gco2_other5 + gco2_comb5
    r5_g_cf4 = gcf4_5

    route_map = {
        "R1_central_storage": (pulse_series(r1_g_co2), None, None),
        "R2_recirculation": (None, pulse_series(r2_g_va), None),
        "R3_municipal_incineration": (pulse_series(r3_g_co2), None, pulse_series(r3_g_cf4)),
        "R4_hazardous_incineration": (pulse_series(r4_g_co2), None, pulse_series(r4_g_cf4)),
        "R5_perfluoro_company": (pulse_series(r5_g_co2), None, pulse_series(r5_g_cf4)),
    }

    if route not in route_map:
        raise ValueError(f"Unknown route {route}")

    co2_rate, va_rate, cf4_rate = route_map[route]

    emissions = {}
    for i in range(n_samples):
        series: Dict[str, np.ndarray] = {}
        if co2_rate is not None:
            series["CO2"] = co2_rate[i]
        if va_rate is not None:
            # Map to DES/SEV/ISO keys used by climate model
            series[{"des":"DES","sev":"SEV","iso":"ISO"}[a]] = va_rate[i]
        if cf4_rate is not None:
            series["CF4"] = cf4_rate[i]
        emissions[f"sample_{i:05d}"] = LegacyEmissions(t_years=t, emissions_g_per_year=series)

    return emissions


def legacy_route_emissions_continuous(
    litres_total: float,
    agent: str = "des",
    route: str = "R2_recirculation",
    years_emit: float = 50.0,
    n_samples: int = 1000,
    seed: Optional[int] = None,
    years_total: float = 200.0,
    dt_years: float = 0.25,
    cf4_yield_fraction: Optional[float] = None,
) -> Dict[str, LegacyEmissions]:
    """Generate emissions schedules for *continuous* legacy disposal over `years_emit` years.

    This replaces the one-time pulse with a constant-rate disposal/release over `years_emit`,
    after which emissions cease and atmospheric decay continues.

    Parameters
    ----------
    litres_total : float
        Total legacy stock volume to be disposed (litres).
    years_emit : float
        Duration over which the stock is released/destroyed at a constant rate.
    years_total : float
        Total simulation horizon (years).
    dt_years : float
        Time step (years).

    Returns
    -------
    Dict[sample_id, LegacyEmissions]
        Species-resolved emission rates (g/year) on the time grid.
    """
    if litres_total < 0:
        raise ValueError("litres_total must be non-negative")
    if years_emit <= 0:
        raise ValueError("years_emit must be positive")
    if years_total <= 0:
        raise ValueError("years_total must be positive")
    if years_emit > years_total:
        raise ValueError("years_emit must be <= years_total")

    a = agent_key(agent)
    route = route.strip()

    rho_values = {"des": P.rho_d, "sev": P.rho_s, "iso": P.rho_i}
    gfm_values = {"des": P.des_gfm, "sev": P.sev_gfm, "iso": P.iso_gfm}
    carbon_economy_values = {"des": P.des_atecon_C, "sev": P.sev_atecon_C, "iso": P.iso_atecon_C}

    rho = rho_values[a]
    gfm = gfm_values[a]
    carbon_economy = carbon_economy_values[a]

    kg_va_total = litres_total * rho / 1000.0
    g_va_total = kg_va_total * 1000.0

    # Constant disposal rate during emission window
    g_va_per_year = g_va_total / years_emit
    kg_va_per_year = kg_va_total / years_emit

    mol_per_kg = 1000.0 / gfm  # mol/kg

    r = rng(seed)
    van_factor = uniform_samples(r, P.min_van_factor, P.max_van_factor, n_samples)
    network_distance = uniform_samples(r, P.min_network_distance, P.max_network_distance, n_samples)
    mwi_dist = uniform_samples(r, P.min_mwi_dist, P.max_mwi_dist, n_samples)
    pwi_dist = uniform_samples(r, P.min_pwi_dist, P.max_pwi_dist, n_samples)

    demunic = uniform_samples(r, P.min_demunic, P.max_demunic, n_samples)
    depharma = uniform_samples(r, P.min_depharma, P.max_depharma, n_samples)
    deplasma = uniform_samples(r, P.min_deplasma, P.max_deplasma, n_samples)

    munic_other = uniform_samples(r, P.min_munic_other, P.max_munic_other, n_samples)
    pharam_other = uniform_samples(r, P.min_pharam_other, P.max_pharam_other, n_samples)
    plasma_other = uniform_samples(r, P.min_plasma_other, P.max_plasma_other, n_samples)

    # Time grid
    t = np.arange(int(np.round(years_total/dt_years)) + 1) * dt_years
    n_t = len(t)
    active = (t < years_emit).astype(float)  # 1 during emissions, 0 after

    # Route contributions expressed as g/year arrays per sample
    # R1 central storage: transport-only consolidation
    x1 = network_distance
    x2 = van_factor
    r1_perkg = (x1 * x2) / 1000.0  # kgCO2/kgVA
    r1_g_co2_per_year = (r1_perkg * kg_va_per_year) * 1000.0  # gCO2/year during active window

    # R2 recirculation: emit parent VA mass
    r2_g_va_per_year = np.full(n_samples, g_va_per_year)

    MW_CO2 = 44.01
    MW_CF4 = 88.00
    mol_per_year = mol_per_kg * kg_va_per_year  # mol/year during active window

    def combustion_products(de: np.ndarray):
        de_eff = de if cf4_yield_fraction is None else np.full_like(de, cf4_yield_fraction)
        g_co2 = (mol_per_year * (1.0 - de_eff) * carbon_economy * MW_CO2)
        g_cf4 = (mol_per_year * de_eff * carbon_economy * MW_CF4)
        return g_co2, g_cf4

    def transport_co2(dist: np.ndarray):
        # dist * tonneVA/year * van_factor -> kgCO2/year -> gCO2/year
        tonneVA_per_year = kg_va_per_year / 1000.0
        kgco2 = dist * tonneVA_per_year * van_factor
        return kgco2 * 1000.0

    # R3 municipal incineration
    gco2_comb3, gcf4_3 = combustion_products(demunic)
    gco2_trans3 = transport_co2(mwi_dist)
    gco2_other3 = (kg_va_per_year * munic_other) * 1000.0
    r3_g_co2_per_year = gco2_trans3 + gco2_other3 + gco2_comb3
    r3_g_cf4_per_year = gcf4_3

    # R4 hazardous incineration
    gco2_comb4, gcf4_4 = combustion_products(depharma)
    gco2_trans4 = transport_co2(pwi_dist)
    gco2_other4 = (kg_va_per_year * pharam_other) * 1000.0
    r4_g_co2_per_year = gco2_trans4 + gco2_other4 + gco2_comb4
    r4_g_cf4_per_year = gcf4_4

    # R5 plasma/perfluoro company
    gco2_comb5, gcf4_5 = combustion_products(deplasma)
    gco2_trans5 = transport_co2(network_distance)
    gco2_other5 = (kg_va_per_year * plasma_other) * 1000.0
    r5_g_co2_per_year = gco2_trans5 + gco2_other5 + gco2_comb5
    r5_g_cf4_per_year = gcf4_5

    def build_rate_series(per_year: np.ndarray) -> np.ndarray:
        # per_year shape (n_samples,)
        return (per_year[:, None] * active[None, :])

    route_map = {
        "R1_central_storage": ("CO2", build_rate_series(r1_g_co2_per_year), None),
        "R2_recirculation": ({ "des":"DES","sev":"SEV","iso":"ISO"}[a], build_rate_series(r2_g_va_per_year), None),
        "R3_municipal_incineration": ("CO2", build_rate_series(r3_g_co2_per_year), build_rate_series(r3_g_cf4_per_year)),
        "R4_hazardous_incineration": ("CO2", build_rate_series(r4_g_co2_per_year), build_rate_series(r4_g_cf4_per_year)),
        "R5_perfluoro_company": ("CO2", build_rate_series(r5_g_co2_per_year), build_rate_series(r5_g_cf4_per_year)),
    }

    if route not in route_map:
        raise ValueError(f"Unknown route {route}")

    spec_key, main_rate, cf4_rate = route_map[route]

    emissions = {}
    for i in range(n_samples):
        series: Dict[str, np.ndarray] = {}
        series[spec_key] = main_rate[i]
        if cf4_rate is not None:
            series["CF4"] = cf4_rate[i]
        emissions[f"sample_{i:05d}"] = LegacyEmissions(t_years=t, emissions_g_per_year=series)

    return emissions
