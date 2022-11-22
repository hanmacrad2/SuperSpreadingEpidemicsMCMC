#RUN MCMC CREDIBLE INTERVALS

#OUTPUT
OUTPUT_FOLDER = 
  
#PARAMS
num_days = 110
alpha_vec = seq(from = 0.7, to = 3, length = 20)

#MCMC FOR EACH VALUE
par(mfrow = c(4,5))
for(i in seq_along(alpha_vec)){
  print(paste0('i: ', i))
  
  #1. SIMULATE DATA
  print(alpha_vec[i]);
  r0 = alpha_vec[i] + 0.5
  data_sse = SIMULATION_SSE(alpha_vec[i])
  plot.ts(data_sse, main = print(paste0(round(r0, 2))))
}

#MCMC SSE POISSON COMPOUND
#START MCMC
Rprof(tmp <- tempfile())
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(canadaX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sse_output$time_elap = time_elap
Rprof()
summaryRprof(tmp)