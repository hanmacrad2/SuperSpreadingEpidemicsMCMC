#Apply model evidence - model comparion basic

#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"

#***********************
#* METHOD 1: HARMONIC MEAN
#**********************

#***********************
# 1. RUN BASE MCMC
#**********************
list_log_ev_base = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER,FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                                                SSIB = FALSE, SSIC = FALSE)) 
#Plot RESULTS
PLOT_MODEL_EV_RESULTS(list_log_ev_base, model_type = 'Baseline')
mean(list_log_ev_base)
sd(list_log_ev_base)

#***********************
# 2. RUN SSEB MCMC
#**********************
list_log_ev_sseb = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER) 
mean(list_log_ev_sseb)
sd(list_log_ev_sseb)

#Plot
PLOT_MODEL_EV_RESULTS(list_log_ev_sseb)

#***********************
#* METHOD 2: IMPORTANCE SAMPLING -- MODEL EVIDENCE RESULTS
#**********************

#1. BASE
model_type = 'BASE'; run_number = 1
FOLDERX = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
lme_base_is = readRDS(file = paste0(FOLDERX, 'phat_ests_base_', run_number))

#Plot
mean(lme_base_is)
sd(lme_base_is)
PLOT_MODEL_EV_RESULTS(lme_base_is)

#2. SSEB
#INSPECT MODEL EVIDENCE RESULTS
model_type = 'SSEB'; run_number = 1
FOLDERX = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
lme_sseb_is = readRDS(file = paste0(FOLDERX, 'phat_ests_sseb_', run_number))

#Plot
mean(lme_sseb_is)
sd(lme_sseb_is)
PLOT_MODEL_EV_RESULTS(lme_sseb_is)








#INSPECT SINGLE SET OF RESULTS
#Inspect MCMC results
model_type = 'SSEB'; run_number = 1
RESULTS_FOLDER = paste0(CURRENT_OUTPUT_FOLDER, model_type, '/run_', run_number, '/')

#MCMC RESULTS
i = 1
mcmc_sseb_1 = readRDS(file = paste0(RESULTS_FOLDER, 'mcmc_sseb_', i))

lme1 = LOG_MODEL_EVIDENCE(mcmc_sseb_1$log_like_vec)
lme1
