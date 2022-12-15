#PLOT SSEC MCMC RESULTS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/PLOT_SSEC_GRID.R')
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/ssec"

#PARAMETERS 
seedX = 1 #SAME AS MCMC RUN
n_mcmc = 100000 #SAME AS MCMC RUN

#CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)

#LOAD DATA & RESULTS
epi_data_ssec = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb_', seedX, '.rds' ))
mcmc_ssec_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', seedX, '.rds' ))

#PLOT MCMC RESULTS
loglike4 =  LOG_LIKE_SSEC(x, lambda_vec, c(1.2, 0.16)) #-127.2522
df_results_ssec = PLOT_SSEC_MCMC_GRID(epi_data_ssec, mcmc_ssec_output, n_mcmc, seedX, loglike4)
