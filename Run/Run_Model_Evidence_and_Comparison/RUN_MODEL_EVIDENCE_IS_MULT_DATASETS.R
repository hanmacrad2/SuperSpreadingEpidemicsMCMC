#**********************************************
#GET MODEL EVIDENCE VIA IMPORTANCE SAMPLING (2018 Paper)
#***********************
library(SuperSpreadingEpidemicsMCMC)
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/PART_2/SSEB_DATA/"
#RUN NUM
run = 1; n_repeats = 30
OUTER_FOLDER = paste0(OUTER_FOLDER, 'run_', run, '/')

#***********************
# 1. DATA 
#**********************
#BASE_DATA = FALSE; SSEB_DATA = TRUE; NZ_DATA = FALSE

#DATASETS; MATRIX
file_name = 'matrix_data_sseb.rds'
matrix_data = readRDS(file = paste0(OUTER_FOLDER, file_name))

#***************************
# 2. LOAD MCMC & GET MULTIPLE PHAT (log)
#***************************

#*************************
#2a. BASELINE
#*************************
list_is_log_ev_base = LOAD_MCMC_GET_P_HAT_II(matrix_data_sseb, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                          FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                              SSIB = FALSE, SSNB = FALSE))
mean(list_is_log_ev_base)
sd(list_is_log_ev_base)

#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_base)

#*************************
#2c. SSNB
#*************************
list_is_log_ev_ssnb = LOAD_MCMC_GET_P_HAT_II(matrix_data_sseb, OUTER_FOLDER, run = run, n_repeats = n_repeats,
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
list_is_log_ev_sseb = LOAD_MCMC_GET_P_HAT_II(matrix_data_sseb, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                          FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                              SSIB = FALSE, SSIC = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_sseb)
sd(list_is_log_ev_sseb)


