#Apply Model Comparison
library(coda)
library(SuperSpreadingEpidemicsMCMC)

#1. RUN SSI MODEL
mcmc_ssi_output = SSI_MCMC_ADAPTIVE(sim_data_canadaX1)

#2. RUN SSE MODEL
mcmc_sse_output = SSE_MCMC_ADAPTIVE(canadaX)

#3. GET HARMONIC MEANS
bf = get_bayes_factor_harmonic_means(mcmc_ssi_output$log_like_vec, mcmc_sse_output$log_like_vec)
print(paste0('bayes factor = ', bf))

#****************************
#PLOT OUTPUT

#1. SSI
df_ssi = PLOT_SS_MCMC_GRID(sim_data_canadaX1, mcmc_ssi_output) 

#2. SSE
df_sse = PLOT_SS_MCMC_GRID(canadaX, mcmc_sse_output,
                              mcmc_specs = list(n_mcmc = 500000, burn_in_size = 0.05,
                                                mod_start_points = list(m1 = 0.8, m2 = 0.05, m3 = 10),
                                                mod_par_names = c('alpha', 'beta', 'gamma'),
                                                model_type = 'SSE', seed_count = 1, thinning_factor = 10),
                              priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                                 c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                              FLAGS_LIST = list(SSI = FALSE, BURN_IN = TRUE, THIN = TRUE,
                                                DATA_AUG = TRUE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                                PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
