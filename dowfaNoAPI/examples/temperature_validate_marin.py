"""Reproduce Marin-style temperature curves (continuous emissions and 50y stop) for CO2 and desflurane,
and add a destruction+CF4 sensitivity sweep using the legacy pulse framework.

This script is intended for Supplementary validation/illustration.
"""
import numpy as np
import matplotlib.pyplot as plt
from dowfa_src.climate import temperature_from_emissions, DEFAULT_GASES

def build_rate(t, annual_g, stop_year=None):
    r = np.zeros_like(t)
    if stop_year is None:
        r[:] = annual_g
    else:
        r[t < stop_year] = annual_g
    return r

def main():
    years_total = 100.0
    dt = 0.25
    t = np.arange(int(round(years_total/dt))+1) * dt

    # Baselines (as per your Marin reproduction)
    E_CO2 = 92.7e6 * 1e6  # g/yr
    E_DES = 13200 * 1e6   # g/yr

    # Panel A: continuous
    out_co2_A = temperature_from_emissions(t, {"CO2": build_rate(t, E_CO2)})
    out_des_A = temperature_from_emissions(t, {"DES": build_rate(t, E_DES)})

    # Panel B: 50y stop
    out_co2_B = temperature_from_emissions(t, {"CO2": build_rate(t, E_CO2, 50.0)})
    out_des_B = temperature_from_emissions(t, {"DES": build_rate(t, E_DES, 50.0)})

    # Desflurane decay -> CF4 formation sweep: treat CF4 as emitted proportional to DES emissions (illustrative)
    alphas = [0.001, 0.01, 0.1, 0.5]
    fig, axs = plt.subplots(1,2, figsize=(13,5))

    axs[0].plot(t, out_co2_A["T"], lw=2, label="CO₂ 92.7 Mt/y")
    axs[0].plot(t, out_des_A["T"], lw=2, label="Des 13,200 t/y")
    for a in alphas:
        out = temperature_from_emissions(t, {"DES": build_rate(t, E_DES), "CF4": build_rate(t, E_DES*a)})
        axs[0].plot(t, out["T"], lw=2, ls="--", label=f"Des + CF₄ (α={a:g})")
    axs[0].set_title("Panel A: continuous emissions")
    axs[0].set_xlabel("Year"); axs[0].set_ylabel("ΔT (°C)")
    axs[0].grid(True, alpha=0.3); axs[0].legend(fontsize=9)

    axs[1].plot(t, out_co2_B["T"], lw=2, label="CO₂ 50y then stop")
    axs[1].plot(t, out_des_B["T"], lw=2, label="Des 50y then stop")
    for a in alphas:
        out = temperature_from_emissions(t, {"DES": build_rate(t, E_DES, 50.0), "CF4": build_rate(t, E_DES*a, 50.0)})
        axs[1].plot(t, out["T"], lw=2, ls="--", label=f"Des + CF₄ (α={a:g})")
    axs[1].set_title("Panel B: stop after 50 years")
    axs[1].set_xlabel("Year"); axs[1].set_ylabel("ΔT (°C)")
    axs[1].grid(True, alpha=0.3); axs[1].legend(fontsize=9)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
