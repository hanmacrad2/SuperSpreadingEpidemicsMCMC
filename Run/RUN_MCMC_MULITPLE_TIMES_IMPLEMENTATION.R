#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
CURRENT_OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#***********************
# 1. BASE DATA (RUN AUTOMATICALLY)
#**********************

LOC_BASE_DATA = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(LOC_BASE_DATA, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
runX = 1

RUN_MULTIPLE_MCMC_SSEC(data_baseI, CURRENT_OUTPUT_FOLDER) 
