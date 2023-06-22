#RUN BASELINE MCMC

#DATA
file_name = 'data_baseline.rds'
data_baseline = readRDS(paste0(DATA_FOLDER, file_name))
plot.ts(data_baseline)


#2.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_baseline = MCMC_INFER_BASELINE(epi_data_sseb, n_mcmc)