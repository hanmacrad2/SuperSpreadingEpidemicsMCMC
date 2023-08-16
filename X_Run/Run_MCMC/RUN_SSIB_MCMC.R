#RUN SSIB MODEL
n_mcmc = 500

start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssib = MCMC_INFER_SSIB(EPI_DATA)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssib$time_elap = time_elap

#PLOT
PLOT_SSB_MCMC_GRID(EPI_DATA, mcmc_ssib, n_mcmc = n_mcmc,
                   FLAGS_MODELS = list(SSEB = FALSE, SSIB = TRUE),
                   sim_vals = list(m1 = 0.8, m2 = 0.1, m3 = 10))

#EDTS
mcmc_ssib$log_like_vec = mcmc_ssib$log_like_vec[1:2400-1]
