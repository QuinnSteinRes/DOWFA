"""Plot temperature trajectories for WVA disposal routes (R6-R11) on the Marin & Kleinberg scale.

We rescale the WVA activity so that the *mean annual DES mass used* during the active emission
window equals 13,200 t/year, matching the Marin comparison case. This scaling is computed
from the R6_release activity schedule and then applied to all routes so that comparisons
reflect the same underlying anaesthetic activity level.

Outputs:
- Panel A: continuous emissions (no stop) for 100 years
- Panel B: emissions for 50 years then stop; simulation shown out to 200 years

Includes baseline curves:
- CO2 continuous: 92.7 Mt/y
- DES continuous: 13,200 t/y
"""

import numpy as np
import matplotlib.pyplot as plt

from dowfa_src import wva_route_emissions
from dowfa_src.climate import temperature_from_emissions

def build_constant_rate(t, annual_g, stop_year=None):
    r = np.zeros_like(t, dtype=float)
    if stop_year is None:
        r[:] = annual_g
    else:
        r[t < stop_year] = annual_g
    return r

def extend_to_grid(series_short, t_full):
    out = np.zeros_like(t_full, dtype=float)
    n = min(len(series_short), len(out))
    out[:n] = series_short[:n]
    return out

def compute_scale_factor_from_release(years_emit, N=300, seed=42, initial_population_m=45.0, des_percent_mac=12.0):
    """Scaling so mean annual DES mass during [0, years_emit) equals 13,200 t/y.

    We compute this from R6_release (which contains DES emissions), then reuse the same
    scaling factor for destruction routes where DES may not appear as an emitted species.
    """
    target_des_g_per_year = 13200.0 * 1e6  # t/y -> g/y
    sims = wva_route_emissions(
        route="R6_release",
        total_years=years_emit,
        initial_population_m=initial_population_m,
        des_percent_mac=des_percent_mac,
        n_samples=N,
        seed=seed,
    )
    any_sample = next(iter(sims.values()))
    des_rate = any_sample.emissions_g_per_year["DES"]
    t_short = any_sample.t_years
    active = (t_short < years_emit)
    mean_des = float(np.mean(des_rate[active]))
    if mean_des <= 0:
        raise ValueError("Mean DES emission is non-positive; cannot scale.")
    return target_des_g_per_year / mean_des

def scaled_temperature_for_route(route, years_emit, continuous, scale, N=300, seed=42,
                                 years_total=200.0, dt=1.0,
                                 initial_population_m=45.0, des_percent_mac=12.0):
    sims = wva_route_emissions(
        route=route,
        total_years=years_emit,
        initial_population_m=initial_population_m,
        des_percent_mac=des_percent_mac,
        n_samples=N,
        seed=seed,
    )

    t_full = np.arange(0.0, years_total + 1e-12, dt)
    Ts = []

    for em in sims.values():
        emissions_full = {}
        for spec, rate_short in em.emissions_g_per_year.items():
            if continuous:
                rate_full = np.zeros_like(t_full, dtype=float)
                n = min(len(rate_short), len(rate_full))
                rate_full[:n] = rate_short[:n]
                if n < len(rate_full):
                    # hold last value constant
                    rate_full[n:] = rate_short[-1]
            else:
                rate_full = extend_to_grid(rate_short, t_full)

            emissions_full[spec] = rate_full * scale

        out = temperature_from_emissions(t_full, emissions_full)
        Ts.append(out["T"])

    Ts = np.vstack(Ts)
    return t_full, np.median(Ts, axis=0), np.percentile(Ts, 5, axis=0), np.percentile(Ts, 95, axis=0)

def main():
    routes = [
        "R6_release",
        "R7_AC+MSWI",
        "R8_AC+PWI",
        "R9_AC+WM",
        "R10_AC+onsite_inc",
        "R11_in_situ_plasma",
    ]
    years_emit = 50
    N = 300
    seed = 42
    initial_population_m = 45.0
    des_percent_mac = 12.0

    scale = compute_scale_factor_from_release(
        years_emit=years_emit, N=N, seed=seed, initial_population_m=initial_population_m, des_percent_mac=des_percent_mac
    )

    dt = 1.0
    years_total_A = 100.0
    years_total_B = 200.0
    tA = np.arange(0.0, years_total_A + 1e-12, dt)
    tB = np.arange(0.0, years_total_B + 1e-12, dt)

    E_CO2 = 92.7e6 * 1e6   # g/yr
    E_DES = 13200.0 * 1e6  # g/yr

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: continuous emissions baselines
    out_co2_A = temperature_from_emissions(tA, {"CO2": build_constant_rate(tA, E_CO2)})
    out_des_A = temperature_from_emissions(tA, {"DES": build_constant_rate(tA, E_DES)})

    axs[0].plot(tA, out_co2_A["T"], lw=2, label="CO2 92.7 Mt/y (Marin)")
    axs[0].plot(tA, out_des_A["T"], lw=2, label="DES 13,200 t/y (Marin)")

    for route in routes:
        t, med, lo, hi = scaled_temperature_for_route(
            route=route, years_emit=years_emit, continuous=True, scale=scale,
            N=N, seed=seed, years_total=years_total_A, dt=dt,
            initial_population_m=initial_population_m, des_percent_mac=des_percent_mac
        )
        axs[0].plot(t, med, lw=2, label=f"{route} (scaled)")
        axs[0].fill_between(t, lo, hi, alpha=0.15)

    axs[0].set_title("Panel A: continuous emissions (Marin scale)")
    axs[0].set_xlabel("Year")
    axs[0].set_ylabel("Delta T (degC)")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(fontsize=8)

    # Panel B: stop after 50 years baselines
    out_co2_B = temperature_from_emissions(tB, {"CO2": build_constant_rate(tB, E_CO2, stop_year=50.0)})
    out_des_B = temperature_from_emissions(tB, {"DES": build_constant_rate(tB, E_DES, stop_year=50.0)})

    axs[1].plot(tB, out_co2_B["T"], lw=2, label="CO2 50y then stop (Marin)")
    axs[1].plot(tB, out_des_B["T"], lw=2, label="DES 50y then stop (Marin)")

    for route in routes:
        t, med, lo, hi = scaled_temperature_for_route(
            route=route, years_emit=years_emit, continuous=False, scale=scale,
            N=N, seed=seed, years_total=years_total_B, dt=dt,
            initial_population_m=initial_population_m, des_percent_mac=des_percent_mac
        )
        axs[1].plot(t, med, lw=2, label=f"{route} (scaled)")
        axs[1].fill_between(t, lo, hi, alpha=0.15)

    axs[1].set_title("Panel B: emissions 0-50y then stop")
    axs[1].set_xlabel("Year")
    axs[1].set_ylabel("Delta T (degC)")
    axs[1].grid(True, alpha=0.3)
    axs[1].legend(fontsize=8)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
