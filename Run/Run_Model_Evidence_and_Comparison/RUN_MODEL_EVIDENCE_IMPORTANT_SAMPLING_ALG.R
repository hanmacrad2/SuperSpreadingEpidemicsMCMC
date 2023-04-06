#IMPLEMENT MODEL EVIDENCE VIA IMPORTANCE SAMPLING 
library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASE_DATA/')

#***********************
# 1. BASE DATA (RUN AUTOMATICALLY)
#**********************
#LOC_BASE_DATA = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(OUTER_FOLDER, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
runX = 1

#***************************
# 2. LOAD MCMC & GET MULTIPLE PHAT (log)
#***************************

#*************************
#2a. BASELINE
#*************************
ests_phat_base = LOAD_MCMC_GET_P_HAT(data_baseI, OUTER_FOLDER, run = 2, n_repeats = 500,
                                     FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                         SSIB = FALSE, SSNB = FALSE))
mean(ests_phat_base)
sd(ests_phat_base)

#PLOT
PLOT_MODEL_EV_RESULTS(ests_phat_base)

#LOAD estimates
run = 2
model_type = 'BASELINE'
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/') #phat_ests_base_', run, '.rds' ) 
ests_phat_base = readRDS(file = paste0(CURRENT_FOLDER, 'phat_ests_base_', run, '.rds' ))

#*************************
#2b. SSEB
#*************************
model_type = 'SSEB'
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/')

ests_phat_sseb = LOAD_MCMC_GET_P_HAT(data_baseI, OUTER_FOLDER,
                                     FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                                                  SSIB = FALSE, SSNB = FALSE))
#PLOT
PLOT_MODEL_EV_RESULTS(ests_phat_sseb)

#*************************
#2c. SSNB
#*************************
ests_phat_ssnb = LOAD_MCMC_GET_P_HAT(data_baseI, OUTER_FOLDER,
                                     FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                          SSIB = FALSE, SSIC = FALSE))
#PLOT
PLOT_MODEL_EV_RESULTS(ests_phat_ssnb)

#Results
mean(ests_phat_ssnb)
sd(ests_phat_ssnb)

