"""
Validation benchmark: Marin & Kleinberg (SI Section IV) pulse experiment.

Reproduces the "pulse" temperature-response comparison:
- Desflurane pulse: 13,200 t at t=0
- CO2 pulse equivalents: 13,200 t * GWP20 and 13,200 t * GWP100

Overlays digitized Marin curves provided as three CSV files:
- marin_black_7020x_co2.csv
- marin_blue_2590x_co2.csv
- marin_red_des.csv

Axis limits follow the figure:
- min_x=min_y=0, max_x=100, max_y=0.00007

Output:
- ../FIGURES/benchmark_marin_pulse.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from dowfa_src.climate import temperature_from_emissions

# ----------------------------
# Plot styling (larger fonts)
# ----------------------------
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 13,
    "axes.linewidth": 1.2,
})

# ----------------------------
# Digitized data (your CSVs)
# ----------------------------
DIGITIZED_BLACK_7020 = "../MartinKleinbergDataset/marin_black_7020x_co2.csv"
DIGITIZED_BLUE_2590  = "../MartinKleinbergDataset/marin_blue_2590x_co2.csv"
DIGITIZED_RED_DES    = "../MartinKleinbergDataset/marin_red_des.csv"


def _load_two_col_csv(path: str):
    """
    Loads a 2-column CSV (year, y) and returns (x, y) as float arrays.
    Robust to headers with leading/trailing spaces.
    """
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding=None)

    # Handle potential whitespace in column names (your files have leading spaces)
    names = [n.strip() for n in data.dtype.names]
    if len(names) != 2:
        raise ValueError(f"Expected exactly 2 columns in {path}, found {names}")

    x = np.array(data[data.dtype.names[0]], dtype=float)
    y = np.array(data[data.dtype.names[1]], dtype=float)
    return x, y


def build_pulse(t, mass_g, at_year=0.0):
    """
    Represent a pulse emission as a one-year impulse at the closest time index.
    Since our convolution is annual (dt=1y), we approximate delta(t) by placing
    all mass in the first bin.
    """
    e = np.zeros_like(t, dtype=float)
    idx = int(np.argmin(np.abs(t - at_year)))
    e[idx] = mass_g
    return e


def summarize_curve(t, T, name):
    i_peak = int(np.argmax(T))
    print(f"{name:>22s}: peak ΔT = {T[i_peak]:.6g} °C at year {t[i_peak]:.0f}")


def main():
    # ----------------------------
    # Time grid
    # ----------------------------
    years_total = 100
    dt = 1.0
    t = np.arange(0.0, years_total + 1e-12, dt)

    # ----------------------------
    # Marin benchmark inputs
    # ----------------------------
    # Pulse mass (tons -> g)
    m_des_tons = 13200.0
    m_des_g = m_des_tons * 1e6  # 1 ton = 1e6 g

    # Marin figure labels (as used in the paper figure)
    GWP20_DES = 7020.0
    GWP100_DES = 2590.0  # NOTE: Marin figure shows 2590x (not 2560x)

    m_co2_gwp20_g = m_des_g * GWP20_DES
    m_co2_gwp100_g = m_des_g * GWP100_DES

    # ----------------------------
    # Emission arrays (pulse at t=0)
    # ----------------------------
    E_des = build_pulse(t, m_des_g, at_year=0.0)
    E_co2_20 = build_pulse(t, m_co2_gwp20_g, at_year=0.0)
    E_co2_100 = build_pulse(t, m_co2_gwp100_g, at_year=0.0)

    # ----------------------------
    # Temperature responses (model)
    # ----------------------------
    out_des = temperature_from_emissions(t, {"DES": E_des})
    out_co2_20 = temperature_from_emissions(t, {"CO2": E_co2_20})
    out_co2_100 = temperature_from_emissions(t, {"CO2": E_co2_100})

    summarize_curve(t, out_des["T"], "DES pulse (model)")
    summarize_curve(t, out_co2_20["T"], "CO2 7020x (model)")
    summarize_curve(t, out_co2_100["T"], "CO2 2590x (model)")

    # ----------------------------
    # Load digitized Marin curves
    # ----------------------------
    xb, yb = _load_two_col_csv(DIGITIZED_BLACK_7020)
    xu, yu = _load_two_col_csv(DIGITIZED_BLUE_2590)
    xr, yr = _load_two_col_csv(DIGITIZED_RED_DES)

    # ----------------------------
    # Plot (Marin = lines, Model = dots)
    # ----------------------------
    fig = plt.figure(figsize=(9.0, 6.0))

    # --- Marin digitised curves (LINES) ---
    plt.plot(
        xb, yb,
        color="black", linewidth=2.8, linestyle="-",
        label=r"7020$\times$13,200 tons CO$_2$, Marin & Kleinberg"
    )

    plt.plot(
        xu, yu,
        color="blue", linewidth=2.8, linestyle="-",
        label=r"2590$\times$13,200 tons CO$_2$, Marin & Kleinberg"
    )

    plt.plot(
        xr, yr,
        color="red", linewidth=2.8, linestyle="-",
        label=r"13,200 tons Desflurane, Marin & Kleinberg"
    )


    # --- Your model curves (DOTS) ---
    plt.scatter(
        t, out_co2_20["T"],
        s=28, color="black", marker="o", alpha=0.9,
        label=r"7020$\times$13,200 tons CO$_2$"
    )

    plt.scatter(
        t, out_co2_100["T"],
        s=28, color="blue", marker="o", alpha=0.9,
        label=r"2590$\times$13,200 tons CO$_2$"
    )

    plt.scatter(
        t, out_des["T"],
        s=28, color="red", marker="o", alpha=0.9,
        label=r"13,200 tons Desflurane"
    )


    plt.xlabel("Years")
    plt.ylabel(r"Global Temperature Change ($^\circ$C)")

    plt.xlim(0.0, 100.0)
    plt.ylim(0.0, 7.0e-5)

    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", frameon=True)
    plt.tight_layout()

    #plt.savefig("../FIGURES/benchmark_marin_pulse.png",
    #            dpi=600, bbox_inches="tight")
    plt.savefig("../FIGURES/benchmark_marin_pulse.eps",
                    dpi=600, bbox_inches="tight", format="eps")

if __name__ == "__main__":
    main()