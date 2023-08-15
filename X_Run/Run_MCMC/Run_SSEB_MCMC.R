#RUN SSEB MCMC 
#FOLDERS
library(SuperSpreadingEpidemicsMCMC)                                           

#OUTPUT FOLDER
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/sseb"
seedX = 7

#CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

#PARAMETERS 
n_mcmc = 30000 
alphaX = 0.8; num_days = 50 #110 #Question: 50 days ok? Run 110 now. 1000 runs takes x mins

#1.SIMULATE DATA
set.seed(seedX)
epi_data_sseb = SIMULATE_EPI_SSEB(num_days)
plot.ts(epi_data_sseb)
#SAVE DATA
saveRDS(epi_data_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb_', seedX, '.rds' ))
epi_data_sseb = sim7

#2.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sseb = MCMC_INFER_SSEB(data_sseb, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sseb$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_sseb_output7, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', seedX, '.rds' ))

#SEED
seedX = seedX + 1

#RUN MCMC
n_mcmc = 50000
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sseb2 = MCMC_INFER_SSEB(data_sseb2, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sseb2$time_elap = time_elap

#Plot
PLOT_SSB_MCMC_GRID(data_sseb, mcmc_sseb, n_mcmc = n_mcmc,
                   FLAGS_MODELS = list(SSEB = TRUE, SSIB = FALSE),
                   sim_vals = list(m1 = 0.8, m2 = 0.1, m3 = 10))
#PLOT_SSB_MCMC_GRIDPLOT_SSB_MCMC_GRID(data_sseb, mcmc_sseb, n_mcmc)
