"""Temperature-response model (AGTP-style impulse-response).

This module implements a lightweight temperature model suitable for Supplementary Material:
- CO2 uses a 4-term impulse response function (IRF) from Joos et al. / IPCC style.
- Other gases use single-exponential lifetime decay.
- Radiative forcing is computed via radiative efficiency (RE) converted to W/m^2 per gram.
- Temperature is computed by convolving RF(t) with a 2-box climate response kernel.

Units:
- Time is in years
- Emissions rates are in g/year on the discretised time grid
- Burdens are in g (atmospheric mass anomaly)
- RF is W/m^2
- Temperature response is in degrees C (relative)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np

# --- Atmosphere conversion ---
M_ATM_G = 5.1352e21  # g
MW_AIR = 28.97       # g/mol

def A_from_RE(RE_wm2_per_ppb: float, MW_g_per_mol: float) -> float:
    """Convert radiative efficiency (W/m^2/ppb) to W/m^2 per gram emitted in the atmosphere."""
    return RE_wm2_per_ppb * (MW_AIR / MW_g_per_mol) * (1e9 / M_ATM_G)

# --- Default parameters (can be overridden) ---
@dataclass(frozen=True)
class CO2IRF:
    a0: float = 0.2173
    a1: float = 0.2240
    t1: float = 394.4
    a2: float = 0.2824
    t2: float = 36.54
    a3: float = 0.2763
    t3: float = 4.304

def f_co2(t: np.ndarray, irf: CO2IRF = CO2IRF()) -> np.ndarray:
    t = np.asarray(t)
    return (irf.a0 +
            irf.a1*np.exp(-t/irf.t1) +
            irf.a2*np.exp(-t/irf.t2) +
            irf.a3*np.exp(-t/irf.t3))

def f_exp(t: np.ndarray, tau: float) -> np.ndarray:
    t = np.asarray(t)
    return np.exp(-t/tau)

@dataclass(frozen=True)
class TempKernel:
    c1: float = 0.631
    d1: float = 8.4
    c2: float = 0.429
    d2: float = 409.5

def h_temp(t: np.ndarray, k: TempKernel = TempKernel()) -> np.ndarray:
    t = np.asarray(t)
    return (k.c1/k.d1)*np.exp(-t/k.d1) + (k.c2/k.d2)*np.exp(-t/k.d2)

@dataclass
class GasSpec:
    name: str
    mw_g_per_mol: float
    re_wm2_per_ppb: float
    lifetime_years: Optional[float] = None  # None => CO2 special

    @property
    def A(self) -> float:
        return A_from_RE(self.re_wm2_per_ppb, self.mw_g_per_mol)

# A small default set aligned with the Marin-style reproduction you were using.
DEFAULT_GASES: Dict[str, GasSpec] = {
    # Values consistent with Marin & Kleinberg supplementary material for CO2 and desflurane
    "CO2": GasSpec("CO2", 44.01, 1.33e-5, None),
    "DES": GasSpec("DES", 168.04, 0.464, 13.9),  # desflurane ~13.9–14.1 years
    # Andersen et al. report radiative efficiencies for VAs (W m−2 ppb−1)
    "SEV": GasSpec("SEV", 200.05, 0.351, 1.7),
    "ISO": GasSpec("ISO", 184.49, 0.453, 3.5),
    # IPCC AR6 (Chapter 7) gives methane radiative efficiency for CH4-fossil:
    "CH4": GasSpec("CH4", 16.04, 5.7e-4, 11.8),
    # CF4 radiative efficiency commonly assessed around 0.085 W m−2 ppb−1 (recent measurements & assessments);
    # lifetime is extremely long (~50,000 years).
    "CF4": GasSpec("CF4", 88.00, 0.085, 50000.0),
}


def temperature_from_emissions(
    t: np.ndarray,
    emissions_g_per_year: Dict[str, np.ndarray],
    gases: Dict[str, GasSpec] = DEFAULT_GASES,
    dt_years: Optional[float] = None,
    co2_irf: CO2IRF = CO2IRF(),
    kernel: TempKernel = TempKernel(),
) -> Dict[str, np.ndarray]:
    """Compute RF(t) and DeltaT(t) from time series emissions per gas.

    Parameters
    ----------
    t : array
        Time grid in years (monotonic).
    emissions_g_per_year : dict
        Map gas -> emission rate array (same length as t), in g/year.
    gases : dict
        Gas specification map; must include the keys present in emissions_g_per_year.
    dt_years : float|None
        Time step; inferred from t if None.

    Returns
    -------
    dict with keys: "RF_total", "T", and optionally per-gas RF arrays.
    """
    t = np.asarray(t)
    n = len(t)
    if dt_years is None:
        if n < 2:
            raise ValueError("t must have at least 2 points")
        dt_years = float(t[1]-t[0])

    # Precompute lags
    lags = np.arange(n) * dt_years
    decay_cache = {}
    for gas, spec in gases.items():
        if gas == "CO2":
            decay_cache[gas] = f_co2(lags, co2_irf)
        else:
            if spec.lifetime_years is None:
                raise ValueError(f"Gas {gas} must have lifetime_years unless CO2")
            decay_cache[gas] = f_exp(lags, spec.lifetime_years)

    # Burden-based RF from convolution of past emissions with decay
    RF_gas = {gas: np.zeros(n) for gas in emissions_g_per_year.keys()}
    for gas, rate in emissions_g_per_year.items():
        rate = np.asarray(rate)
        if len(rate) != n:
            raise ValueError(f"Emission series for {gas} length mismatch")
        A = gases[gas].A
        decay = decay_cache[gas]
        for i in range(n):
            dm = rate[:i+1] * dt_years  # g emitted in each step
            RF_gas[gas][i] = A * np.sum(dm[::-1] * decay[:i+1])

    RF_total = np.zeros(n)
    for gas, rf in RF_gas.items():
        RF_total += rf

    # Temperature convolution
    h = h_temp(lags, kernel)
    T = np.zeros(n)
    for i in range(n):
        T[i] = np.sum(RF_total[:i+1][::-1] * h[:i+1]) * dt_years

    out = {"RF_total": RF_total, "T": T}
    for gas, rf in RF_gas.items():
        out[f"RF_{gas}"] = rf
    return out
