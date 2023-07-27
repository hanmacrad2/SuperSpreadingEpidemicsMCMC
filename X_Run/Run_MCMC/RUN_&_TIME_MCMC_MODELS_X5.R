#RUN & TIME EACH MODEL
library(SuperSpreadingEpidemicsMCMC)
n_mcmc = 30000

#1. BASELINE MODEL
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_baseline = MCMC_INFER_BASELINE(EPI_DATA, n_mcmc)

#2.SSE (NB)
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_sse = MCMC_INFER_SSNB(EPI_DATA, n_mcmc)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
print(time_elap)

#3. SSEB
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_sseb = MCMC_INFER_SSEB(EPI_DATA, n_mcmc)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
print(time_elap)

#4. SSIB
print(paste0('start_time:', start_time))

mcmc_ssib = MCMC_INFER_SSIB(EPI_DATA, n_mcmc = n_mcmc)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
print(time_elap)

#5.SSI
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_ssi = MCMC_INFER_SSIR(EPI_DATA, n_mcmc = n_mcmc)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
print(time_elap)