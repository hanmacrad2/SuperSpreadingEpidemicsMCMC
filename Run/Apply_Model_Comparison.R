#Apply Model Comparison
library(coda)
#1. RUN SSI MODEL
mcmc_ssi_output = SSI_MCMC_ADAPTIVE(sim_data_canadaX1)

#2. RUN SSE MODEL
mcmc_sse_output = SSE_MCMC_ADAPTIVE(canadaX)

#3. GET HARMONIC MEANS
bf = get_bayes_factor_harmonic_means(mcmc_ssi_output$log_like_vec, mcmc_sse_output$log_like_vec)
print(paste0('bayes factor = ', bf))

#PLOT OUTPUT
df_ssi = PLOT_SS_MCMC_GRID(sim_data_canadaX1, mcmc_ssi_output) 
