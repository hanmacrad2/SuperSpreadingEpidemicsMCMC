#**********************************************
#GET MODEL EVIDENCE VIA IMPORTANCE SAMPLING (2018 Paper)
#***********************
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/PART_2/SSEB_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21_SUBSET_I/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21/"
run = 1; 
n_repeats = 50

#***********************
# 1. DATA 
#**********************
#BASE_DATA = FALSE; SSEB_DATA = TRUE; NZ_DATA = FALSE

#BASE DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
file_name = 'epi_data_base_1.rds'
data_baseline = readRDS(file = paste0(OUTER_FOLDER, file_name))
plot.ts(data_baseline)

#SSEB DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
file_name = "epi_data_sseb_1.rds"
data_sseb = readRDS(file = paste0(OUTER_FOLDER, file_name))
plot.ts(data_sseb)

#NZ DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21/"
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
data_file_wait_21 = read.csv(paste0(DATA_FOLDER, 'data_waitemata_aug_21.csv'))
data_wait_08_21 = data_file_wait_21$Cases
plot.ts(data_wait_08_21)

#***************************
#***************************

#*************************
#2a. BASELINE
#*************************
list_is_log_ev_base = LOAD_MCMC_GET_P_HAT(data_sseb, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                     FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                         SSIB = FALSE, SSNB = FALSE))
mean(list_is_log_ev_base)
sd(list_is_log_ev_base)

#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_base)

#LOAD estimates
#model_type = 'BASELINE'
#CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/') 
#list_is_log_ev_base = readRDS(file = paste0(CURRENT_FOLDER, 'phat_ests_base_', run, '.rds' ))

#*************************
#2c. SSNB
#*************************
run = 3; n_repeats = 50
list_is_log_ev_ssnb = LOAD_MCMC_GET_P_HAT(data_wait_08_21, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                           FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                               SSIB = FALSE, SSIC = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_ssnb)

#Results
mean(list_is_log_ev_ssnb)
sd(list_is_log_ev_ssnb)

#*************************
#2b. SSEB
#*************************
run = 2
list_is_log_ev_sseb = LOAD_MCMC_GET_P_HAT(data_wait_08_21, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                          FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                              SSIB = FALSE, SSIC = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_sseb)
sd(list_is_log_ev_sseb)

#OUTPUT
model_type = 'SSEB'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/')
file_name = 'phat_ests_sseb_2.rds'
model_ev_sseb =  readRDS(file = paste0(CURRENT_FOLDER, file_name))

#OUTPUT
model_type = 'SSNB'; print(model_type)
run = 3
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/')
file_name = 'phat_ests_ssnb_3.rds'
model_ev_ssnb3 =  readRDS(file = paste0(CURRENT_FOLDER, file_name))

#*****************

#COUNTEIS MANAUKA

#****************
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'CM_08_21_SUB_1/')
print(OUTER_FOLDER)

n_repeats = 34; run =1

DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
filename =  'data_cm_tot_aug_21.csv'

data_file_cm_21 = read.csv(paste0(DATA_FOLDER, filename))
data_file_cm_21 = data_file_cm_21$Cases
plot.ts(data_file_cm_21, ylab = 'Infection count', main = 'CM NZ, August 2021')

#SUBSET
data_file_cm_21_sub1 = data_file_cm_21[5:21]
plot.ts(data_file_cm_21_sub1)

#*************************
#2a. BASELINE
#*************************
list_is_log_ev_baseCM = LOAD_MCMC_GET_P_HAT(data_file_cm_21_sub1, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                            start = 16,
                                          FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                              SSIB = FALSE, SSNB = FALSE))
mean(list_is_log_ev_baseCM, na.rm = TRUE)
sd(list_is_log_ev_baseCM, na.rm = TRUE)

#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_base)

#LOAD estimates
#model_type = 'BASELINE'
#CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/') 
#list_is_log_ev_base = readRDS(file = paste0(CURRENT_FOLDER, 'phat_ests_base_', run, '.rds' ))

#*************************
#2c. SSNB
#*************************
run = 1; 
list_is_log_ev_ssnbCM = LOAD_MCMC_GET_P_HAT(data_file_cm_21_sub1, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                            start = 16,
                                          FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                              SSIB = FALSE, SSIC = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_ssnb)

#Results
mean(list_is_log_ev_ssnbCM)
sd(list_is_log_ev_ssnbCM)

#*************************
#2b. SSEB
#*************************
run = 1
list_is_log_ev_sseb = LOAD_MCMC_GET_P_HAT(data_wait_08_21_sub1, OUTER_FOLDER, run = 1, n_repeats = n_repeats,
                                            start = 1,
                                          FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                              SSIB = FALSE, SSIC = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_ssebCM)
sd(list_is_log_ev_ssebCM)

#*************************
#2b. SSIR
#*************************
run = 4; n_repeats = 10
list_is_log_ev_ssir = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub1, OUTER_FOLDER,
                                            run = run, n_repeats = n_repeats,
                                            FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                SSIB = FALSE, SSIR = TRUE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_ssir)
sd(list_is_log_ev_ssir)

#*************
#LOAD MODEL EVIDENCE ESTIMATES
#**************
run = 1
model_type = 'baseline'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
list_is_log_ev_base = readRDS(file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))

