#RUN SSEB MCMC 
#FOLDERS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R'); source('R/EPI_FUNCTIONS.R')
source('R/SSEC_MCMC.R')

#OUTPUT FOLDER
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/ssec"
seedX = 4

#CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

#PARAMETERS 
n_mcmc = 10000 #50000 
#alphaX = 0.8; num_days = 50 #110 #Question: 50 days ok? Run 110 now. 1000 runs takes x mins

#1.SIMULATE DATA
seedX = seedX + 1
set.seed(seedX)
epi_data_ssec = SIMULATE_EPI_SSEC()
plot.ts(epi_data_ssec)
#SAVE DATA
saveRDS(epi_data_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb_', seedX, '.rds' ))
epi_data_sseb = sim7
#2.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssec_output = MCMC_INFER_SSEC(epi_data_ssec, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssec_output$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_ssec_output, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssec_', seedX, '.rds' ))

#SEED
seedX = seedX + 1