#RUN SSIB MODEL

start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssib = MCMC_INFER_SSIB(EPI_DATA)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssib$time_elap = time_elap