"""
DES-only Marin-symmetry plots with Monte-Carlo uncertainty bands.

Key points:
- Scale factor computed from R6_release so mean DES during years 0–50 equals 13,200 t/y.
- SEV/ISO are excluded (set to zero) to match Marin's DES-only reference.
- For destruction routes, CO2/CF4 by-products are retained (otherwise routes collapse to ~0).
- Monte-Carlo bands are computed from the ensemble of temperature trajectories:
    median = p50, band = p5–p95.
- Panel A: continuous emissions (extend after 50y by holding last-year rate constant), shown to 100y.
- Panel B: emissions for 0–50y then stop, shown to 200y.
- Marin baselines plotted: CO2 (black) and DES (red), unshaded.
- EPS-safe styling: no transparency, solid light uncertainty bands, publication-style formatting.
"""

import numpy as np
import matplotlib.pyplot as plt

from dowfa_src import wva_route_emissions
from dowfa_src.climate import temperature_from_emissions


# ----------------------------
# Publication-style plotting setup
# ----------------------------
plt.rcParams.update({
    "font.size": 15,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "legend.fontsize": 13,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "ps.fonttype": 42,   # better EPS font handling
    "pdf.fonttype": 42
})

KEEP_SPECS = {"DES", "CO2", "CF4"}  # drop SEV/ISO only

P_LO = 5.0
P_HI = 95.0

ROUTE_LABELS = {
    "R6_release": "R6: Direct release",
    "R7_AC+MSWI": "R7: AC + MSWI",
    "R8_AC+PWI": "R8: AC + PWI",
    "R9_AC+WM": "R9: AC + specialist waste",
    "R10_AC+onsite_inc": "R10: AC + on-site incineration",
    "R11_in_situ_plasma": "R11: In-situ plasma",
}

# EPS-safe publication colours
ROUTE_COLORS = {
    "R6_release": "#D55E00",      # orange-red
    "R7_AC+MSWI": "#0072B2",      # blue
    "R8_AC+PWI": "#009E73",       # green
    "R9_AC+WM": "#CC79A7",        # purple
    "R10_AC+onsite_inc": "#E69F00",
    "R11_in_situ_plasma": "#56B4E9",
}

ROUTE_BAND_COLORS = {
    "R6_release": "#F6C7A6",
    "R7_AC+MSWI": "#BFD9EE",
    "R8_AC+PWI": "#BFE3D4",
    "R9_AC+WM": "#E7C6DB",
    "R10_AC+onsite_inc": "#F3D59B",
    "R11_in_situ_plasma": "#CFEAF7",
}


def build_constant_rate(t, annual_g, stop_year=None):
    r = np.zeros_like(t, dtype=float)
    if stop_year is None:
        r[:] = annual_g
    else:
        r[t < stop_year] = annual_g
    return r


def compute_scale_factor(years_emit, N=500, seed=42,
                         initial_population_m=45.0, des_percent_mac=12.0):
    """Scale so mean annual DES over [0, years_emit) equals 13,200 t/y."""
    target_des = 13200.0 * 1e6  # g/year
    sims = wva_route_emissions(
        route="R6_release",
        total_years=years_emit,
        initial_population_m=initial_population_m,
        des_percent_mac=des_percent_mac,
        n_samples=N,
        seed=seed,
    )
    sample = next(iter(sims.values()))
    des_rate = sample.emissions_g_per_year["DES"]
    mean_des = float(np.mean(des_rate[:years_emit]))
    if mean_des <= 0:
        raise ValueError("Non-positive DES mean; cannot scale.")
    return target_des / mean_des


def route_temperature_band(route, years_emit, continuous, scale, years_total,
                           dt=1.0, N=500, seed=42,
                           initial_population_m=45.0, des_percent_mac=12.0):
    """Return (t, median, lo, hi) for a route."""
    sims = wva_route_emissions(
        route=route,
        total_years=years_emit,
        initial_population_m=initial_population_m,
        des_percent_mac=des_percent_mac,
        n_samples=N,
        seed=seed,
    )

    t = np.arange(0.0, years_total + 1e-12, dt)
    Ts = []

    for em in sims.values():
        emissions = {}
        for spec, rate in em.emissions_g_per_year.items():
            if spec not in KEEP_SPECS:
                continue

            full = np.zeros_like(t)
            n = min(len(rate), len(t))
            full[:n] = rate[:n]

            if continuous and n < len(t):
                full[n:] = rate[-1]

            emissions[spec] = full * scale

        out = temperature_from_emissions(t, emissions)
        Ts.append(out["T"])

    Ts = np.vstack(Ts)
    med = np.percentile(Ts, 50.0, axis=0)
    lo = np.percentile(Ts, P_LO, axis=0)
    hi = np.percentile(Ts, P_HI, axis=0)
    return t, med, lo, hi


def main():

    routes = [
        "R6_release",
        "R7_AC+MSWI",
        "R8_AC+PWI",
        # "R9_AC+WM",
        # "R10_AC+onsite_inc",
        # "R11_in_situ_plasma",
    ]

    years_emit = 50
    N = 500
    seed = 42
    initial_population_m = 45.0
    des_percent_mac = 12.0

    scale = compute_scale_factor(
        years_emit, N=N, seed=seed,
        initial_population_m=initial_population_m,
        des_percent_mac=des_percent_mac
    )

    # Marin baselines
    E_CO2 = 92.7e6 * 1e6   # g/yr
    E_DES = 13200.0 * 1e6  # g/yr

    # ========== Panel A: continuous ==========
    tA = np.arange(0.0, 100.0 + 1e-12, 1.0)
    out_co2_A = temperature_from_emissions(tA, {"CO2": build_constant_rate(tA, E_CO2)})
    out_des_A = temperature_from_emissions(tA, {"DES": build_constant_rate(tA, E_DES)})

    plt.figure(figsize=(9, 5.5))

    plt.plot(
        tA, out_co2_A["T"],
        color="black", linewidth=2.4,
        label="CO₂ 92.7 Mt/y (Marin & Kleinberg, 2025)",
        zorder=3
    )
    plt.plot(
        tA, out_des_A["T"],
        color="red", linewidth=2.4,
        label="Desflurane 13,200 t/y (Marin & Kleinberg, 2025)",
        zorder=3
    )

    for r in routes:
        t, med, lo, hi = route_temperature_band(
            r, years_emit, continuous=True, scale=scale, years_total=100.0,
            N=N, seed=seed,
            initial_population_m=initial_population_m,
            des_percent_mac=des_percent_mac
        )

        label_text = ROUTE_LABELS.get(r, r)
        line_color = ROUTE_COLORS.get(r, "0.3")
        band_color = ROUTE_BAND_COLORS.get(r, "0.85")

        plt.fill_between(
            t, lo, hi,
            facecolor=band_color,
            edgecolor="none",
            zorder=1
        )
        plt.plot(
            t, med,
            color=line_color,
            linewidth=2.2,
            label=label_text,
            zorder=2
        )

    plt.xlabel("Year")
    plt.ylabel("ΔT (°C)")

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(
        loc="best",
        frameon=False,
        handlelength=2.8,
        labelspacing=0.4,
        borderpad=0.2
    )
    plt.tight_layout()
    plt.savefig("../FIGURES/contOnlyDeswva.eps", dpi=600, bbox_inches="tight", format="eps")

    # ========== Panel B: 0–50y then stop ==========
    tB = np.arange(0.0, 200.0 + 1e-12, 1.0)
    out_co2_B = temperature_from_emissions(tB, {"CO2": build_constant_rate(tB, E_CO2, stop_year=50.0)})
    out_des_B = temperature_from_emissions(tB, {"DES": build_constant_rate(tB, E_DES, stop_year=50.0)})

    plt.figure(figsize=(9, 5.5))

    plt.plot(
        tB, out_co2_B["T"],
        color="black", linewidth=2.4,
        label="CO₂ 50y then stop (Marin & Kleinberg, 2025)",
        zorder=3
    )
    plt.plot(
        tB, out_des_B["T"],
        color="red", linewidth=2.4,
        label="Desflurane 50y then stop (Marin & Kleinberg, 2025)",
        zorder=3
    )

    for r in routes:
        t, med, lo, hi = route_temperature_band(
            r, years_emit, continuous=False, scale=scale, years_total=200.0,
            N=N, seed=seed,
            initial_population_m=initial_population_m,
            des_percent_mac=des_percent_mac
        )

        label_text = ROUTE_LABELS.get(r, r)
        line_color = ROUTE_COLORS.get(r, "0.3")
        band_color = ROUTE_BAND_COLORS.get(r, "0.85")

        plt.fill_between(
            t, lo, hi,
            facecolor=band_color,
            edgecolor="none",
            zorder=1
        )
        plt.plot(
            t, med,
            color=line_color,
            linewidth=2.2,
            label=label_text,
            zorder=2
        )

    plt.xlabel("Year")
    plt.ylabel("ΔT (°C)")

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(
        loc="best",
        frameon=False,
        handlelength=2.8,
        labelspacing=0.4,
        borderpad=0.2
    )
    plt.tight_layout()
    plt.savefig("../FIGURES/cease50yearsOnlyDeswva.eps", dpi=600, bbox_inches="tight", format="eps")


if __name__ == "__main__":
    main()