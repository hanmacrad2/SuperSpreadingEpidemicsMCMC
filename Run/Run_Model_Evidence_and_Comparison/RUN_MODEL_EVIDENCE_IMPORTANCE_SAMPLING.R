#**********************************************
#GET MODEL EVIDENCE VIA IMPORTANCE SAMPLING (2018 Paper)
#***********************
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/PART_2/SSEB_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21_SUBSET_I/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21_SUBSET_II/"
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

#*************************************************
#*
#* MODEL EVIDENCE ESTIMATES
#* 
#* ***************************************************************************

#*************************
#2a. BASELINE
#*************************
run = 1
list_is_log_ev_base = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub1, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = FALSE))
mean(list_is_log_ev_base, na.rm = TRUE)
sd(list_is_log_ev_base, na.rm = TRUE)

#*************************
#2c. SSNB
#*************************
run = 5; 
list_is_log_ev_ssnb = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub1, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                        SSIB = FALSE, SSIR = FALSE))

#Results
mean(list_is_log_ev_ssnb)
sd(list_is_log_ev_ssnb)

#*************************
#2b. SSEB
#*************************
run = 1
list_is_log_ev_sseb = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub2, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = FALSE))
#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_sseb)
sd(list_is_log_ev_sseb)

#*************************
#2b. SSIR
#*************************
run = 5
list_is_log_ev_ssir = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub1, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = TRUE))

#PLOT
#PLOT_MODEL_EV_RESULTS(list_is_log_ev_sseb)
mean(list_is_log_ev_ssir)
sd(list_is_log_ev_ssir)

#*************************
#2b. SSNB - GAMMA PRIOR ON K
#*************************
run = 'gamma_prior_k'; 
list_is_log_ev_ssnb_gak = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub2, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                       FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                           SSIB = FALSE, SSIR = FALSE))

PLOT_MODEL_EV_RESULTS(list_is_log_ev_ssnb_gak)

#PLOT
data_type = 'Waitemata days 10-20'
par(mfrow = c(1,1))
BOX_PLOT_MODEL_EV(list_vec_results = list(ssnb_k_exp_prior = list_is_log_ev_ssnb,
                                          ssnb_k_ga_prior = list_is_log_ev_ssnb_gak),
                              title = 'SSNB Model k priors compared exp(1) vs Ga(0.001, rt = 0.001). ',
                              data_type = data_type, model = '') 

#*****************************
#*
#LOAD MODEL EVIDENCE ESTIMATES
#*
#*****************************

run = 1
model_type = 'baseline'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
list_is_log_ev_base = readRDS(file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))