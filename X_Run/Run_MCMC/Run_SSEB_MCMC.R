#RUN SSEB MCMC 
n_mcmc = 30000
#RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sseb3 = MCMC_INFER_SSEB(data_sseb3, n_mcmc = n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sseb3$time_elap = time_elap

#SAVE
file1 = 'mcmc_sseb_d2_50k.rds'
filename = 'mcmc_sseb_d1_50k.rds'
mcmc_sseb2 = readRDS(paste0(OUTER_FOLDER, filename))

saveRDS(mcmc_sseb, paste0(OUTER_FOLDER, file1))

#Plot
PLOT_SSB_MCMC_GRID(data_sseb3, mcmc_sseb3, n_mcmc = n_mcmc,
                   FLAGS_MODELS = list(SSEB = TRUE, SSIB = FALSE),
                   sim_vals = list(m1 = 0.7, m2 = 0.05, m3 = 10))   
                   #sim_vals = list(m1 = 0.9, m2 = 0.05, m3 = 8))                                      

#SAVE
file1 = 'mcmc_sseb_09_05_8'
saveRDS(mcmc_sseb, paste0(OUTER_FOLDER, file1))


#RUN ORIG CODE
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sseb11 = MCMC_INFER_SSEB(data_sseb)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)

#Plot
PLOT_SSB_MCMC_GRID(data_sseb, mcmc_sseb11, n_mcmc = n_mcmc,
                   FLAGS_MODELS = list(SSEB = TRUE, SSIB = FALSE),
                   sim_vals = list(m1 = 0.9, m2 = 0.05, m3 = 8))                                      

#SAVE
file1 = 'mcmc_sseb_09_05_8_orig_p2'
saveRDS(mcmc_sseb11, paste0(OUTER_FOLDER, file1))
mcmc_sseb10 = readRDS(paste0(OUTER_FOLDER, file1))

mcmc_sseb11$alpha_vec = mcmc_sseb11$alpha_vec[1:length(mcmc_sseb11$alpha_vec)-1]

mcmc_sseb11$beta_vec = mcmc_sseb11$beta_vec[1:length(mcmc_sseb11$beta_vec)-1]

mcmc_sseb11$gamma_vec = mcmc_sseb11$gamma_vec[1:length(mcmc_sseb11$gamma_vec)-1]

mcmc_sseb11$r0_vec = mcmc_sseb11$r0_vec[1:length(mcmc_sseb11$r0_vec)-1]

mcmc_sseb11$log_like_vec = mcmc_sseb11$log_like_vec[1:length(mcmc_sseb11$log_like_vec)-1]




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


