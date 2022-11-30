#RUN MCMC CREDIBLE INTERVALS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R')
source('R/SSE_MCMC_POISSON_COMPOUND.R')

#PARAMETERS
alphaX = 1.2
data_sse = SIMULATION_SSE(alphaX)
plot.ts(data_sse)
