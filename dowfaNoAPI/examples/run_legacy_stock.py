from dowfa_src import simulate_legacy_stock

res = simulate_legacy_stock(litres=250.0, agent="des", gwp="GWP_100", n_samples=10000, seed=42)

print("Totals (tonnes CO2e):", res.totals_tonnes_co2e)
print("Per-kg (kgCO2e/kg VA):", res.per_kg_kgco2e_per_kgVA)
print("Std dev per-kg:", res.stdev_per_kg)
