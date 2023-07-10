#RUN BASELINE MCMC
n_mcmc = 30000

#DATA
file_name = 'data_baseline.rds'
data_baseline = readRDS(paste0(DATA_FOLDER, file_name))
plot.ts(data_baseline)

#2.RUN MCMC - exp prior
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_baseline = MCMC_INFER_BASELINE(EPI_DATA, n_mcmc)

mcmc_baseline$r0_vec = head(mcmc_baseline$r0_vec, -1)

mcmc_baseline$log_like_vec = head(mcmc_baseline$log_like_vec, -1)
                            
#PLOT
PLOT_BASELINE_R0_MCMC(EPI_DATA, mcmc_baseline)

#2.RUN MCMC - exp prior
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_baseline_ga = MCMC_INFER_BASELINE(data_baseline, n_mcmc,
                                    FLAGS_LIST = list(ADAPTIVE = TRUE, 
                                                      PRIOR_EXP = TRUE, PRIOR_GAMMA = TRUE,
                                                      THIN = FALSE, BURN_IN = TRUE))

PLOT_BASELINE_R0_MCMC(data_baseline, mcmc_baseline_ga)
