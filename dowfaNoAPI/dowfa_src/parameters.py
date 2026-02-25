"""Parameter definitions for DOWFA calculation library.

"""

# GWP values
des_20 = 6975
sev_20 = 604
iso_20 = 1925
cf4_20 = 4400

des_100 = 2560
sev_100 = 168
iso_100 = 538
cf4_100 = 6500

des_500 = 731
sev_500 = 48
iso_500 = 154
cf4_500 = 10000

# Liquid densities (kg/m^3)
rho_d = 1440
rho_s = 1505
rho_i = 1450

# Molecular weights (g/mol)
des_gfm = 168.04
sev_gfm = 200.05
iso_gfm = 184.49

# Atomic economy (C atoms per molecule, used in your original combustion proxy)
des_atecon_C = 3
sev_atecon_C = 4
iso_atecon_C = 3

# Destruction efficiency ranges (fraction forming CF4 proxy in your original model)
min_demunic = 0.17
max_demunic = 0.40

min_depharma = 0.03
max_depharma = 0.10

min_deplasma = 0.0001
max_deplasma = 0.001

# Other emissions factors (kgCO2e per kg waste / process proxy)
min_munic_other = 0.7
max_munic_other = 1.7

min_pharam_other = 0.851
max_pharam_other = 1.05

min_plasma_other = 0.040
max_plasma_other = 0.089

# MAC / surgery assumptions
MAC_des = 6.7
MAC_sev = 2.2
MAC_iso = 1.2

min_sur_100k = 10332
max_sur_100k = 10590

min_growth_rate = 0.0034
max_growth_rate = 0.0036

# 2017 baseline
sur_pro_base = 4.54  # million surgeries
ktco2e_des_base = 89.25
ktco2e_sev_base = 10.5
ktco2e_iso_base = 5.25

# Activated carbon / GAC48
cococnut_to_ac_ratio = 6.7
ac_co2e_per_kg = 1.15

indo_to_sing_port = 979
sing_port_to_fac = 27
sing_to_eng_port = 15205

mol_per_kg_des = 2.63
mol_per_kg_sev = 2.95
mol_per_kg_iso = 1.53

# Transport factors
min_van_factor = 0.13926
max_van_factor = 0.29416

min_hgv_factor = 0.01816
max_hgv_factor = 0.12465

min_ship_factor = 0.00057
max_ship_factor = 0.00833

# Distances (km)
min_network_distance = 2500
max_network_distance = 3500

min_mwi_dist = 24
max_mwi_dist = 48

min_pwi_dist = 93
max_pwi_dist = 193
