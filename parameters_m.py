# parameters.py



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

rho_d = 1440
rho_s = 1505
rho_i = 1450

des_gfm = 168.04
sev_gfm = 200.05
iso_gfm = 184.49

des_atecon_C = 3
sev_atecon_C = 4  
iso_atecon_C = 3  # Atomic Economy for Carbon

#Emmisions Data--------------------------------------

#demunic = 0.285
min_demunic = 0.17
max_demunic = 0.40

#depharma = 0.065
min_depharma = 0.03
max_depharma = 0.10

#deplasma = 0.00055
min_deplasma = 0.0001
max_deplasma = 0.001

#CO2e / kg waste 
#munic_other = 1.2 #Assuming dillution is not required and just dosed 
min_munic_other = 0.7
max_munic_other = 1.7

#pharam_other = 0.96
min_pharam_other = 0.851
max_pharam_other = 1.05

#plasma_other = 0.389 #Estimated from, electrical use 
min_plasma_other = 0.040
max_plasma_other = 0.089

#Other------------------------------------------
MAC_des = 6.7
MAC_sev = 2.2
MAC_iso = 1.2

#sur_100k = 10487 # Intermediate and restrictive procedures per 100k ave UK average 2009-2014
min_sur_100k = 10332
max_sur_100k = 10590


#growth_rate = 0.0035 # growthrate 
min_growth_rate = 0.0034
max_growth_rate = 0.0036


#2017 Base-----------------------------------------
mol_2017_des = 209
mol_2017_sev = 404
mol_2017_iso = 56

sur_pro_base = 4.54 # M surgeries in 2017
sur_emm_base = 105 # ktonneCO2e from anaesthetic gas use in 2017

per_des_sur_base = 0.85 # share of des 2017
per_sev_sur_base = 0.10 # share of des 2017
per_iso_sur_base = 0.05 # share of des 2017

mol_des_base = 209 #Cumulative mol required for surgery in 2017
mol_sev_base = 404 #Cumulative mol required for surgery in 2017
mol_iso_base = 56 #Cumulative mol required for surgery in 2017

ktco2e_des_base = 89.25
ktco2e_sev_base = 10.5
ktco2e_iso_base = 5.25


#AC Params

cococnut_to_ac_ratio = 6.7
ac_co2e_per_kg = 1.15

indo_to_sing_port = 979 #km
sing_port_to_fac = 27 #km
sing_to_eng_port = 15205 #km

#GAC48------------------------------------------------

mol_per_kg_des = 2.63
mol_per_kg_sev = 2.95
mol_per_kg_iso = 1.53

#Transportation Factors ------------------------------

#van_factor = 0.2169
min_van_factor = 0.13926
max_van_factor = 0.29416

#hgv_factor = 0.0714
min_hgv_factor = 0.01816
max_hgv_factor = 0.12465

#ship_factor = 0.00445
min_ship_factor = 0.00057
max_ship_factor = 0.00833


#Transportation Distance--------------------------------

#network_distance = 3000 # km
min_network_distance = 2500
max_network_distance = 3500

#mwi_dist = 36
min_mwi_dist = 24
max_mwi_dist = 48

#pwi_dist = 143
min_pwi_dist = 93
max_pwi_dist = 193
  
  

  


