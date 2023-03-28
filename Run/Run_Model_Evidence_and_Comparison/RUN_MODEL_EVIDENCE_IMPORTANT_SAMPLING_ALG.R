#IMPLEMENT MODEL EVIDENCE VIA IMPORTANCE SAMPLING 
library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#***********************
# 1. BASE DATA (RUN AUTOMATICALLY)
#**********************
LOC_BASE_DATA = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(LOC_BASE_DATA, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
runX = 1

#***************************
# 2. LOAD MCMC & GET MULTIPLE PHAT (log)
#***************************
#SSEB
OUTPUT_FOLDER = paste0(LOC_BASE_DATA, 'SSEB/')

ests_phat_sseb = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                                                  SSIB = FALSE, SSNB = FALSE))
#PLOT
PLOT_MODEL_EV_RESULTS(ests_phat_sseb)
#SAVE
#saveRDS(ests_phat_sseb, file = paste0(OUTPUT_FOLDER, '/run_', runX, '/ests_phat_sseb.rds'))

#BASELINE
OUTPUT_FOLDER = paste0(LOC_BASE_DATA, 'BASE/')
ests_phat_base2 = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                       SSIB = FALSE, SSNB = FALSE))
mean(ests_phat_base)
sd(ests_phat_base)

#PLOT
PLOT_MODEL_EV_RESULTS(ests_phat_base)

#SAVE
saveRDS(ests_phat_base, file = paste0(OUTPUT_FOLDER, '/run_',
                                      runX, '/ests_phat_base_100_r.rds'))

ests_phat_base1 = ests_phat_base

#OUTPUT_FOLDER = paste0(LOC_BASE_DATA, 'SSIB/')
#ests_phat_ssib = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
#                                     FLAGS_LIST = list(BASE = FALSE, SSEB = FALSE,
#                                                       SSIB = TRUE, SSNB = FALSE))

#LOAD MCMC
i = 10
mcmc_samples = readRDS(file = paste0(OUTPUT_FOLDER, 'run_',
                                     runX, '/mcmc_sseb_', i ))
