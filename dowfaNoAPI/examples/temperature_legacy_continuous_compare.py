"""Plot legacy disposal scenarios: continuous emissions for 50 years then stop.

Produces ΔT(t) curves for R2 release vs destruction routes, including CF4 sensitivity override.
"""
import numpy as np
import matplotlib.pyplot as plt
from dowfa_src import legacy_route_emissions_continuous
from dowfa_src.climate import temperature_from_emissions

def run_route(route, litres_total=250.0, years_emit=50.0, N=500, seed=42, cf4=None):
    sims = legacy_route_emissions_continuous(
        litres_total=litres_total,
        agent="des",
        route=route,
        years_emit=years_emit,
        n_samples=N,
        seed=seed,
        years_total=200.0,
        dt_years=0.25,
        cf4_yield_fraction=cf4,
    )
    # summarize median
    t = None
    Ts = []
    for em in sims.values():
        t = em.t_years
        out = temperature_from_emissions(t, em.emissions_g_per_year)
        Ts.append(out["T"])
    Ts = np.vstack(Ts)
    return t, np.median(Ts, axis=0)

def main():
    routes = ["R2_recirculation", "R3_municipal_incineration", "R4_hazardous_incineration", "R5_perfluoro_company"]
    fig, ax = plt.subplots(figsize=(8,5))
    for route in routes:
        t, T = run_route(route)
        ax.plot(t, T, lw=2, label=route)

    # CF4 sensitivity for plasma route
    for cf4 in [0.001, 0.01, 0.1, 0.5]:
        t, T = run_route("R5_perfluoro_company", cf4=cf4)
        ax.plot(t, T, lw=2, ls="--", label=f"R5 plasma + CF4 y={cf4:g}")

    ax.set_xlabel("Year")
    ax.set_ylabel("ΔT (°C)")
    ax.set_title("Legacy disposal: continuous 50y then stop (median over Monte Carlo)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
