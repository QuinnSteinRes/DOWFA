"""Temperature modelling for WVA routes (R6–R11), with emissions continuing for 50 years then ceasing.

We:
1) build yearly emissions time series for each route (0..50 years)
2) extend them with zeros to a longer horizon (default 500 years) so we can evaluate long-term decay
3) compute ΔT(t) using the impulse-response temperature model.

Outputs quantiles (5/50/95%) for ΔT at 20y, 100y, 500y.
"""
import numpy as np
from dowfa_src import wva_route_emissions
from dowfa_src.climate import temperature_from_emissions

def summarize(samples, q=(5,50,95)):
    return np.percentile(samples, q)

def extend_to_horizon(t_short, series_short, years_total=500.0, dt=1.0):
    """Pad to a uniform grid up to years_total with zeros after last value."""
    t_full = np.arange(0.0, years_total + 1e-12, dt)
    out = np.zeros_like(t_full)
    n = min(len(series_short), len(t_full))
    out[:n] = series_short[:n]
    return t_full, out

def main():
    N = 1000
    seed = 42
    years_emit = 50  # continuous for 50 years then stop
    years_total = 500.0
    dt = 1.0

    routes = [
        "R6_release",
        "R7_AC+MSWI",
        "R8_AC+PWI",
        "R9_AC+WM",
        "R10_AC+onsite_inc",
        "R11_in_situ_plasma",
    ]

    for route in routes:
        sims = wva_route_emissions(
            route=route,
            total_years=years_emit,
            initial_population_m=45.0,
            des_percent_mac=12.0,
            n_samples=N,
            seed=seed,
        )

        T20=[]; T100=[]; T500=[]
        for em in sims.values():
            t_short = em.t_years
            emissions_full = {}
            for spec, rate in em.emissions_g_per_year.items():
                t_full, rate_full = extend_to_horizon(t_short, rate, years_total=years_total, dt=dt)
                emissions_full[spec] = rate_full

            out = temperature_from_emissions(t_full, emissions_full)
            T = out["T"]
            T20.append(T[np.argmin(np.abs(t_full-20.0))])
            T100.append(T[np.argmin(np.abs(t_full-100.0))])
            T500.append(T[np.argmin(np.abs(t_full-500.0))])

        print(f"\n{route} (emit 0–{years_emit}y then stop):")
        print("  ΔT@20y  (5/50/95%):", summarize(np.array(T20)))
        print("  ΔT@100y (5/50/95%):", summarize(np.array(T100)))
        print("  ΔT@500y (5/50/95%):", summarize(np.array(T500)))

if __name__ == "__main__":
    main()
