#RUN SSIB MODE
n_mcmc = 100000
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssib_d3 = MCMC_INFER_SSIB(data_ssib, n_mcmc = n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssib_d3$time_elap = time_elap

#PLOT
#pdf('mcmc_grid_ssib', 10, 10)
sim_vals = list(m1 = a, m2 = b, m3 = c)
PLOT_SSB_MCMC_GRID(data_ssib, mcmc_ssib_d3, n_mcmc = n_mcmc,
                   FLAGS_MODELS = list(SSEB = FALSE, SSIB = TRUE),
                   sim_vals = sim_vals)
dev.off()
#SAVE
file = 'mcmc_ssib_50_1_3_ga_8_1_150k.rds'
file_name = 'mcmc_ssib_d1.rds'
saveRDS(mcmc_ssib_d2, paste0(OUTER_FOLDER, file_name))

#EDTS
#mcmc_ssib2$log_like_vec = mcmc_ssib2$log_like_vec[1:2399-1]
