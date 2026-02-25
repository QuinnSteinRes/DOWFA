from dowfa_src import simulate_wva

res = simulate_wva(total_years=10, initial_population_m=56.0, des_percent_mac=5.0, n_samples=5000, seed=1)

print("Totals (ktonnes CO2e):")
for k,v in res.totals_ktonnes_co2e.items():
    print(f"  {k}: {v:.2f}")

print("\nAnnual (ktonnes CO2e/yr):")
for k,v in res.annual_ktonnes_co2e.items():
    print(f"  {k}: {v:.2f}")
