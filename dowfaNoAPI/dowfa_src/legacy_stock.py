"""Legacy stock CO2e calculator 

Monte Carlo simulations

"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Any
import numpy as np

from . import parameters as P
from .utils import rng, uniform_samples, agent_key, gwp_horizon_key

@dataclass(frozen=True)
class LegacyStockResult:
    agent: str
    gwp_horizon: str
    litres_input: float
    n_samples: int
    totals_tonnes_co2e: Dict[str, float]
    per_kg_kgco2e_per_kgVA: Dict[str, float]
    stdev_per_kg: Dict[str, float]

def simulate_legacy_stock(
    litres: float,
    agent: str = "des",
    gwp: str | int = "GWP_100",
    n_samples: int = 10_000,
    seed: Optional[int] = None,
) -> LegacyStockResult:
    """Simulate CO2e for legacy stock disposal routes.

    Parameters
    ----------
    litres : float
        Volume of liquid VA stock (L).
    agent : str
        'des', 'sev', or 'iso'.
    gwp : str|int
        20, 100, 500 or 'GWP_20' etc.
    n_samples : int
        Monte Carlo sample count.
    seed : int|None
        RNG seed.

    Returns
    -------
    LegacyStockResult
    """
    if litres < 0:
        raise ValueError("litres must be non-negative")

    a = agent_key(agent)
    g = gwp_horizon_key(gwp)

    gwp_categories = {
        "GWP_20": {"des": P.des_20, "sev": P.sev_20, "iso": P.iso_20},
        "GWP_100": {"des": P.des_100, "sev": P.sev_100, "iso": P.iso_100},
        "GWP_500": {"des": P.des_500, "sev": P.sev_500, "iso": P.iso_500},
    }
    rho_values = {"des": P.rho_d, "sev": P.rho_s, "iso": P.rho_i}
    gfm_values = {"des": P.des_gfm, "sev": P.sev_gfm, "iso": P.iso_gfm}
    carbon_economy_values = {"des": P.des_atecon_C, "sev": P.sev_atecon_C, "iso": P.iso_atecon_C}
    cf4_gwp_values = {"GWP_20": P.cf4_20, "GWP_100": P.cf4_100, "GWP_500": P.cf4_500}

    rho = rho_values[a]
    gfm = gfm_values[a]
    carbon_economy = carbon_economy_values[a]
    gwp_va = gwp_categories[g][a]
    cf4_gwp = cf4_gwp_values[g]

    # Original conversion (litres -> kg): kg = litres * (rho kg/m3) / 1000
    kg_va = litres * rho / 1000.0

    # mol per kg (original code: mol_val = 1000/gfm) where gfm is g/mol
    mol_per_kg = 1000.0 / gfm  # mol/kg

    r = rng(seed)

    demunic = uniform_samples(r, P.min_demunic, P.max_demunic, n_samples)
    depharma = uniform_samples(r, P.min_depharma, P.max_depharma, n_samples)
    deplasma = uniform_samples(r, P.min_deplasma, P.max_deplasma, n_samples)

    munic_other = uniform_samples(r, P.min_munic_other, P.max_munic_other, n_samples)
    pharam_other = uniform_samples(r, P.min_pharam_other, P.max_pharam_other, n_samples)
    plasma_other = uniform_samples(r, P.min_plasma_other, P.max_plasma_other, n_samples)

    van_factor = uniform_samples(r, P.min_van_factor, P.max_van_factor, n_samples)
    mwi_dist = uniform_samples(r, P.min_mwi_dist, P.max_mwi_dist, n_samples)
    pwi_dist = uniform_samples(r, P.min_pwi_dist, P.max_pwi_dist, n_samples)
    network_distance = uniform_samples(r, P.min_network_distance, P.max_network_distance, n_samples)

    # R1: Central storage (transport-only consolidation)
    # As per manuscript: R1 = x1 * x2, where
    #   x1 = cumulative collection distance [2500, 3500] km
    #   x2 = van emission factor [0.13962, 0.29416] kgCO2e / tonneVA / km
    # Convert to per-kg VA by dividing by 1000 (tonne -> kg).
    x1 = network_distance
    x2 = van_factor
    r1_perkg = (x1 * x2) / 1000.0  # kgCO2e / kgVA
    r1_t = (r1_perkg * kg_va) / 1000.0  # tonnes CO2e for given kg_va

    # Combustion proxy from your original model:
    # combust_value = CO2 from (1 - de) + CF4 from (de) scaled by cf4_gwp
    # carbon_economy * 44.1 and 88.0043 are molar masses used in your original.
    def combust_value(de: np.ndarray) -> np.ndarray:
        return (((mol_per_kg * (1.0 - de)) * carbon_economy * 44.1) / 1000.0) +                (((mol_per_kg * de) * carbon_economy * 88.0043) / 1000.0) * cf4_gwp

    combust_munic = combust_value(demunic)
    combust_pharma = combust_value(depharma)
    combust_plasma = combust_value(deplasma)
# Route totals (tonnes CO2e)
    # Release: kg_va * gwp_va -> kgCO2e, /1000 -> tonnes
    release_t = (kg_va * gwp_va) / 1000.0

    # Incineration routes include: transport + combustion proxy + other emissions proxy
    mswi_t = ((mwi_dist * (kg_va/1000.0) * van_factor) + (kg_va * combust_munic) + (kg_va * munic_other)) / 1000.0
    pwi_t  = ((pwi_dist * (kg_va/1000.0) * van_factor) + (kg_va * combust_pharma) + (kg_va * pharam_other)) / 1000.0
    plasma_t = ((network_distance * (kg_va/1000.0) * van_factor) + (kg_va * combust_plasma) + (kg_va * plasma_other)) / 1000.0

    # Per-kg metrics (kgCO2e/kg VA)
    release_perkg = (release_t * 1000.0) / kg_va if kg_va > 0 else 0.0
    mswi_perkg = (mswi_t * 1000.0) / kg_va if kg_va > 0 else 0.0
    pwi_perkg = (pwi_t * 1000.0) / kg_va if kg_va > 0 else 0.0
    plasma_perkg = (plasma_t * 1000.0) / kg_va if kg_va > 0 else 0.0

    totals = {
        "R1_central_storage": float(np.mean(r1_t)),
        "R2_recirculation": float(np.mean(release_t) if np.ndim(release_t) else release_t),
        "R3_municipal_incineration": float(np.mean(mswi_t)),
        "R4_hazardous_incineration": float(np.mean(pwi_t)),
        "R5_perfluoro_company": float(np.mean(plasma_t)),
    }

    perkg = {
        "R1_central_storage": float(np.mean(r1_perkg)),
        "R2_recirculation": float(np.mean(release_perkg) if np.ndim(release_perkg) else release_perkg),
        "R3_municipal_incineration": float(np.mean(mswi_perkg)),
        "R4_hazardous_incineration": float(np.mean(pwi_perkg)),
        "R5_perfluoro_company": float(np.mean(plasma_perkg)),
    }

    stdev = {
        "R1_central_storage": float(np.std(r1_perkg)),
        "R3_municipal_incineration": float(np.std(mswi_perkg)),
        "R4_hazardous_incineration": float(np.std(pwi_perkg)),
        "R5_perfluoro_company": float(np.std(plasma_perkg)),
    }

    return LegacyStockResult(
        agent=a,
        gwp_horizon=g,
        litres_input=float(litres),
        n_samples=int(n_samples),
        totals_tonnes_co2e=totals,
        per_kg_kgco2e_per_kgVA=perkg,
        stdev_per_kg=stdev,
    )
