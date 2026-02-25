from dowfa_src import simulate_legacy_stock, simulate_wva

def pm(mean, std, fmt="{:.2f}"):
    return f"{fmt.format(mean)} ± {fmt.format(std)}"

def main():
    N = 1000
    SEED = 42

    legacy = simulate_legacy_stock(litres=1.0, agent="des", gwp="GWP_100", n_samples=N, seed=SEED)

    print("===== TABLE 2: Legacy stock (kgCO2e/kg Desflurane) =====")
    order = [
        "R1_central_storage",
        "R2_recirculation",
        "R3_municipal_incineration",
        "R4_hazardous_incineration",
        "R5_perfluoro_company",
    ]
    for route in order:
        mean = legacy.per_kg_kgco2e_per_kgVA[route]
        std = legacy.stdev_per_kg.get(route, 0.0)
        print(f"{route:28s}  {pm(mean, std)}")

    wva = simulate_wva(total_years=3, initial_population_m=45.0, des_percent_mac=12.0, n_samples=N, seed=SEED)

    print("\n===== TABLE 2: WVA routes (kgCO2e/kg WVA) =====")
    for route, means in wva.perkg_kgco2e_per_kgVA.items():
        stds = wva.stdev_perkg_kgco2e_per_kgVA[route]
        des = pm(means["des"], stds["des"], fmt="{:.1f}")
        sev = pm(means["sev"], stds["sev"], fmt="{:.1f}")
        iso = pm(means["iso"], stds["iso"], fmt="{:.1f}")
        print(f"{route:22s}  Des={des:14s}  Sev={sev:14s}  Iso={iso:14s}")

if __name__ == "__main__":
    main()
