#**************************************
#*MODEL EVIDENCE VIA HARMONIC MEAN
#*************************************
#library(SuperSpreadingEpidemicsMCMC)

#FOLDER - MCMC RESULTS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21/"

#PARAMS
run = 1; n_repeats = 100

#***********************
# 1. RUN BASELINE MCMC
#**********************
list_hm_log_ev_base_nz = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER, run = run,  n_repeats = n_repeats,
                                                FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                                    SSNB = FALSE, SSIB = FALSE, SSIC = FALSE))
#STATS
mean(list_hm_log_ev_base_nz)
sd(list_hm_log_ev_base_nz)

#PLOT RESULTS
PLOT_MODEL_EV_RESULTS(list_hm_log_ev_base, model_type = 'Baseline')

#SAVE
saveRDS(list_hm_log_ev_base, file = paste0(OUTPUT_FOLDER, '/run_', run, '/list_hm_log_ev_base.rds'))

#***********************
# 2. RUN SSEB MCMC
#**********************
n_repeats = 90
list_hm_log_ev_sseb = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                                    SSNB = FALSE, SSIB = FALSE, SSIC = FALSE))
mean(list_hm_log_ev_sseb)
sd(list_hm_log_ev_sseb)

#Plot
PLOT_MODEL_EV_RESULTS(list_hm_log_ev_sseb) 

#***********************
# 3. RUN SSNB MCMC
#**********************
list_hm_log_ev_ssnb_nz = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE,
                                                                    SSNB = TRUE, SSIB = FALSE, SSIC = FALSE))
mean(list_hm_log_ev_ssnb_nz)
sd(list_hm_log_ev_ssnb_nz)

#Plot
PLOT_MODEL_EV_RESULTS(list_hm_log_ev_ssnb_nz)

#***********************
# 4. RUN SSIB MCMC 
#**********************


#***********************
# BAYES FACTORS
#**********************
log_bf_hm = list_log_ev_base - list_log_ev_sseb

#PLOT
PLOT_BAYES_FACTORS(log_bf_hm)
mean(log_bf_hm)

