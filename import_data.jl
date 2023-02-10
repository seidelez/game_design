#assign
T=48
T_half=24
P_pv_max=80
dt=1
n_runs=365
rho_efficiency = 1
rho_production = 1
rev=Float64[]
R_SE=Float64[]
R_GEN=Float64[]
R_GEN_tot=Float64[]
C_L=Float64[]

#load data
xf = XLSX.readxlsx("foglio_2.xlsx")
S1=xf["dati"]
prezzi_pun_long=S1["H2:J8761"]
prezzi_magg_tut_long=S1["E2:G8761"]
#scegli anno di simulazione 1=2020, 2=2021, 3=2022
year = 2
p_pun_long=[prezzi_pun_long[1:8760,year];prezzi_pun_long[1:T_half,year]]
p_mt_long=[prezzi_magg_tut_long[1:8760,year];prezzi_magg_tut_long[1:T_half,year]]

lambda_pun = p_pun_long[1:T]
lambda_tarid = p_mt_long[1:T]

PV_ava_long_a=S1["B2:B8761"]
PV_ava_long_a .*= rho_production 
PV_ava_long=[PV_ava_long_a;PV_ava_long_a[1:T_half]]
P_load_long_a=S1["A2:A8761"]
P_load_long=[P_load_long_a;P_load_long_a[1:T_half]]
p_prem=S1["K2"]
p_are=S1["L2"]

T = 24
lambda_pun = p_pun_long[1:T]
lambda_tarif = p_mt_long[1:T]
lambda_prem = p_prem + p_are

#flag_arb=1 per PUN, 0 per maggior tutela
flag_arb=1

