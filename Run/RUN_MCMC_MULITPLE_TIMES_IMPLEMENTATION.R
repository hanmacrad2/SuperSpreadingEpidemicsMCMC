#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
CURRENT_OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#***********************
# 1. EPIDEMIC DATA (RUN AUTOMATICALLY)
#**********************

LOC_BASE_DATA = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(LOC_BASE_DATA, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
run_number = 1

#***********************
# 2. RUN MCMC
#**********************

RUN_MULTIPLE_MCMC_SSEC(data_baseI, CURRENT_OUTPUT_FOLDER) 

#Inspect MCMC results
model_type = 'SSEC'; i = 100
RESULTS_FOLDER = paste0(CURRENT_OUTPUT_FOLDER, '/', model_type, '/run_', run_number, '/')
mcmc100 = readRDS(file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
ssec100 = mcmc_ssec_output$ssec_params_matrix
plot.ts(ssec100)
