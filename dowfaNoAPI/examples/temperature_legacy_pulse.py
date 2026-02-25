"""Compute temperature trajectories for each legacy stock route (continuous disposal then stop).

Legacy stock is disposed/released continuously over `years_emit` years, after which emissions cease.
Monte Carlo parameter ranges match Table 2.

Outputs summary stats (median + 5-95%) for ΔT at 20y, 100y, 500y.
"""
import numpy as np
from dowfa_src import legacy_route_emissions_continuous
from dowfa_src.climate import temperature_from_emissions

def summarize(samples, q=(5,50,95)):
    return np.percentile(samples, q)

def main():
    litres_total = 250.0
    N = 2000
    seed = 42

    routes = [
        "R1_central_storage",
        "R2_recirculation",
        "R3_municipal_incineration",
        "R4_hazardous_incineration",
        "R5_perfluoro_company",
    ]

    years_total = 500.0
    dt = 0.25

    for route in routes:
        sims = legacy_route_emissions_continuous(
            litres_total=litres_total,
            agent="des",
            route=route,
            n_samples=N,
            seed=seed,
            years_total=years_total,
            dt_years=dt,
        )
        t = None
        T20=[]; T100=[]; T500=[]
        for em in sims.values():
            t = em.t_years
            out = temperature_from_emissions(t, em.emissions_g_per_year)
            T = out["T"]
            # index by nearest time
            T20.append(T[np.argmin(np.abs(t-20.0))])
            T100.append(T[np.argmin(np.abs(t-100.0))])
            T500.append(T[np.argmin(np.abs(t-500.0))])
        T20 = np.array(T20); T100=np.array(T100); T500=np.array(T500)
        print(f"\n{route} (continuous over 50y {litres_total} L desflurane):")
        print("  ΔT@20y  (5/50/95%):", summarize(T20))
        print("  ΔT@100y (5/50/95%):", summarize(T100))
        print("  ΔT@500y (5/50/95%):", summarize(T500))

if __name__ == "__main__":
    main()
