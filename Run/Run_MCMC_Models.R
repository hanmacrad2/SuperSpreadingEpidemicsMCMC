#RUN MODELS

#MCMC SSE POISSON COMPOUND

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(canadaX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sse_output$time_elap = time_elap

#PLOT
df_sse = PLOT_SS_MCMC_GRID(canadaX, mcmc_sse_output)

