#**********************************************
#GET MODEL EVIDENCE VIA IMPORTANCE SAMPLING (2018 Paper)
#***********************
library(SuperSpreadingEpidemicsMCMC)
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
run = 2; n_repeats = 50

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

#SSEB DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21/"
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
data_file_wait_21 = read.csv(paste0(DATA_FOLDER, 'data_waitemata_aug_21.csv'))
data_wait_08_21 = data_file_wait_21$Cases
plot.ts(data_wait_08_21)

#***************************
# 2. LOAD MCMC & GET MULTIPLE PHAT (log)
#***************************

#*************************
#2a. BASELINE
#*************************
list_is_log_ev_base = LOAD_MCMC_GET_P_HAT(data_wait_08_21, OUTER_FOLDER, run = run, n_repeats = n_repeats,
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
run = 1
list_is_log_ev_sseb = LOAD_MCMC_GET_P_HAT(data_wait_08_21, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                          FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                              SSIB = FALSE, SSIC = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_sseb)
sd(list_is_log_ev_sseb)


