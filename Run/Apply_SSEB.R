#RUN MCMC 
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R')
source('R/SSEB_MCMC.R')
source('R/PLOT_SS_GRID.R')
seedX = 1

#PARAMETERS
seedX = seedX + 1
set.seed(seedX)
alphaX = 0.8; num_days = 50
data_sse2 = SIMULATION_SSE(alphaX)
plot.ts(data_sse2)
#dataII = data_sse

#2. MCMC
n_mcmc = 100000 #500000 #19.84 seconds :D 
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(data_sse, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sse_output$time_elap = time_elap
#saveRDS(mcmc_sse_output)

#Plot
df_sse = PLOT_SS_MCMC_GRID(data_sse, mcmc_sse_output, n_mcmc,
                           mcmc_specs = list(model_type = 'SSE',
                                             mod_par_names = c('alpha', 'beta', 'gamma'),
                                             seed_count = 1, burn_in_pc = 0.05, thinning_factor = 10),
                           priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                              c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                           FLAGS_LIST = list(SSI = FALSE, BURN_IN = FALSE, THIN = TRUE,
                                             DATA_AUG = FALSE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                             PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))


#2. MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_outputI = SSE_POI_MCMC_ADAPTIVE(data_sse, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sse_outputI$time_elap = time_elap
#saveRDS(mcmc_sse_output)


#Plot
par(mfrow = c(3,3))
#1. SIMULATE DATA
set.seed(i); print(alpha_vec[i]);
data_sse = SIMULATION_SSE(alpha_vec[i])
#SAVE DATA
saveRDS(data_sse)
