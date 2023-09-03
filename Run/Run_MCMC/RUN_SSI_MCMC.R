#RUN SSIB MODE
n_mcmc = 50000
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_3 = MCMC_INFER_SSI(EPI_DATA, n_mcmc = n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_3$time_elap = time_elap

#PLOT
#pdf('mcmc_grid_ssib', 10, 10)
sim_vals = list(m1 = R0, m2 = k)
PLOT_SSI_MCMC_GRID(EPI_DATA, mcmc_ssi_3, EPI_DATA, 1, 
                   0, n_mcmc, #log_like_sim
                   sim_vals = sim_vals)


dev.off()
#SAVE
file = 'mcmc_ssib_50_1_3_ga_8_1_150k.rds'
file_name = 'mcmc_ssib_d1.rds'
saveRDS(mcmc_ssib_d2, paste0(OUTER_FOLDER, file_name))
