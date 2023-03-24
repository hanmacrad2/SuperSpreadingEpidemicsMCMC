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
OUTPUT_FOLDER = paste0(LOC_BASE_DATA, 'SSEB/')
ests_phat_sseb = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                                                  SSIB = FALSE, SSNB = FALSE))


OUTPUT_FOLDER = paste0(LOC_BASE_DATA, 'BASE/')
ests_phat_base = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                       SSIB = FALSE, SSNB = FALSE))

#OUTPUT_FOLDER = paste0(LOC_BASE_DATA, 'SSIB/')
#ests_phat_ssib = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
#                                     FLAGS_LIST = list(BASE = FALSE, SSEB = FALSE,
#                                                       SSIB = TRUE, SSNB = FALSE))

#SAVE
#saveRDS(ests_phat_sseb, file = paste0(OUTPUT_FOLDER, '/run_', runX, '/ests_phat_sseb.rds'))