#RUN SSIB MODE
n_mcmc = 150000
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssib_150k = MCMC_INFER_SSIB(EPI_DATA, n_mcmc = n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssib_150k$time_elap = time_elap

#PLOT
pdf('mcmc_grid_ssib', 10, 10)
PLOT_SSB_MCMC_GRID(EPI_DATA, mcmc_ssib, n_mcmc = n_mcmc,
                   FLAGS_MODELS = list(SSEB = FALSE, SSIB = TRUE),
                   sim_vals = list(m1 = 0.8, m2 = 0.05, m3 = 10))
dev.off()
#SAVE
file1 = 'mcmc_ssib_50_1_3_ga_8_1_100k.rds'
saveRDS(mcmc_ssib, paste0(OUTER_FOLDER, file1))

#EDTS
#mcmc_ssib2$log_like_vec = mcmc_ssib2$log_like_vec[1:2399-1]
