# RUN SSI SIMULATION

#PLOT SIMULATIONS
RESULTS_FOLDER = '~/GitHub/RUN_MODELS/PLOT_SIMULATIONS/DATA/'

#1. SSI MODEL
FLAG_MODEL = GET_FLAGS_MODELS(SSI = TRUE)
MODEL_COLOR = MODEL_COLORS[3]

#DATA
DATA_FOLDER = '~/GitHub/RUN_MODELS/PLOT_SIMULATIONS/DATA/SSI/data/'


#****************
#PLOT R0
FLAGS_PARAM = GET_PARAM(r0 = TRUE)
list_params_ssi_r0 = list(r0 = c(0.9, 1.5, 2.5),
                          k = rep(0.1, 3))
PLOT_SIMULATION_LIST(list_ssi_r0, list_params_ssi_r0, 
                     FLAG_MODEL, FLAGS_PARAM, 
                     MODEL_COLOR, RESULTS_FOLDER) 

#PLOT k
FLAGS_PARAM = GET_PARAM(k = TRUE)
list_params_ssi_k = list(r0 = rep(2.0, 3),
                         k = c(0.05, 0.2, 0.8))
PLOT_SIMULATION_LIST(list_ssi_k, list_params_ssi_k, 
                     FLAG_MODEL, FLAGS_PARAM, 
                     MODEL_COLOR, RESULTS_FOLDER) 


#*************************
#2.SIMULATE DATA
#*************************
#PARAMS
list_params = list(k = rep(0.05, 9), r0 = rep(1.2, 9))

#1. R0, k = 0.1 (R0: 0.9, 1.5, 2.5 )
list_params_ssi_r0 = list(r0 = c(0.9, 1.5, 2.5),
                          k = rep(0.1, 3))

#i. r0_param = 0.9; k = 0.1
r0_param = 0.9; k = 0.1

SSI_DATA_r01 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r01)

#r0_param = 1.5; k = 0.1
r0_param = 1.5; k = 0.1
SSI_DATA_r04 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r04)

SSI_DATA_r05 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r05)

SSI_DATA_r06 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r06)

#r0_param = 2.5; k = 0.1
r0_param = 2.5; k = 0.1
SSI_DATA_r07 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r07)

SSI_DATA_r08 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r08)

SSI_DATA_r09 = SIMULATE_EPI_SSI(r0 = r0_param,  k = k)$epidemic_data
plot.ts(SSI_DATA_r09)

list_epi_data_ssi_r0 = list(SSI_DATA_r02, SSI_DATA_r04, SSI_DATA_r08,
                            SSI_DATA_r01, SSI_DATA_r06, SSI_DATA_r07, 
                            SSI_DATA_r03, SSI_DATA_r05, SSI_DATA_r09)

#SAVE
filename = 'list_ssi_sim_data_r0.rds'
DATA_FOLDER = '~/Project_Folder/Results_Thesis/2_SIMULATIONS/SSI/data/'
create_folder(DATA_FOLDER)
saveRDS(list_epi_data_ssi_r0, file = paste0(DATA_FOLDER, filename))

#PLOT
FLAGS_PARAM = GET_PARAM(r0 = TRUE)
PLOT_SIMULATION_LIST(list_epi_data_ssi_r0, list_params_ssi_r0, 
                     FLAG_MODEL, FLAGS_PARAM, 
                     MODEL_COLOR, RESULTS_FOLDER) 

#*********************
#* 2. k
#*********************
#1. k (k = 0.05, 0.2, 0.8)
list_params_ssi_k = list(r0 = rep(2.0, 3),
                         k = c(0.05, 0.2, 0.8))
# k = c(rep(0.05, 3), rep(0.2, 3), rep(0.8, 3)))

#i. k = 0.05
k_param = 0.05; r0_param = 2.0

SSI_DATA_k1 = SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k1)

SSI_DATA_k2 = SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k2)

SSI_DATA_k3 = SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k3)

#ii. k = 0.2
k_param = 0.2; r0_param = 2.0

SSI_DATA_k4 = SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k4)

SSI_DATA_k5 = SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k5)

SSI_DATA_k6 = SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k6)

#ii. k = 0.8
k_param = 0.8; r0_param = 2.0

SSI_DATA_k7= SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k7)

SSI_DATA_k8= SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k8)

SSI_DATA_k9= SIMULATE_EPI_SSI(r0 = r0_param, k = k_param)$epidemic_data
plot.ts(SSI_DATA_k9)

list_epi_data_ssi_k = list(SSI_DATA_k1, SSI_DATA_k4, SSI_DATA_k7,
                           SSI_DATA_k2, SSI_DATA_k6, SSI_DATA_k9, 
                           SSI_DATA_k3, SSI_DATA_k5, SSI_DATA_k8)

#PLOT
FLAGS_PARAM = GET_PARAM(k = TRUE)
PLOT_SIMULATION_LIST(list_epi_data_ssi_k, list_params_ssi_k, 
                     FLAG_MODEL, FLAGS_PARAM, 
                     MODEL_COLOR, RESULTS_FOLDER) 

#SAVE
filename = 'list_ssi_sim_data_k.rds' 
DATA_FOLDER = '~/Project_Folder/Results_Thesis/2_SIMULATIONS/SSI/data/'
create_folder(DATA_FOLDER)
saveRDS(list_epi_data_ssi_k, file = paste0(DATA_FOLDER, filename))
