#RUN MCMC 
#FOLDERS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R')
source('R/SSEB_MCMC.R')
source('R/PLOT_SS_GRID.R')

#PARAMETERS 
seedX = 1
n_mcmc = 100000 
alphaX = 0.8; num_days = 50 #Question: 50 days ok? Run 110 now

#1. SIMULATE DATA
seedX = seedX + 1
set.seed(seedX)
data_sseb = SIMULATION_SSE(alphaX)
plot.ts(data_sseb)

#2. RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(data_sse, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sse_output$time_elap = time_elap

#SAVE
saveRDS(mcmc_sse_output)