#PLOT SSEB MCMC RESULTS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/PLOT_SS_GRID.R')
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/sseb"

#PARAMETERS (BEST WAY TO LOAD AND SAVE - PROBABLY A FILE)
seedX = 1 #SAME AS MCMC RUN
n_mcmc = 100000 #SAME AS MCMC RUN

#CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)

#LOAD DATA & RESULTS
epi_data_sseb = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb_', seedX, '.rds' ))
mcmc_sseb_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', seedX, '.rds' ))

#PLOT MCMC RESULTS
df_results_sseb = PLOT_SS_MCMC_GRID(epi_data_sseb, mcmc_sseb_output, n_mcmc,
                           mcmc_specs = list(model_type = 'SSEB',
                                             mod_par_names = c('alpha', 'beta', 'gamma'),
                                             seed_count = seedX, burn_in_pc = 0.05, thinning_factor = 10),
                           priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                              c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                           FLAGS_LIST = list(SSI = FALSE, BURN_IN = FALSE, THIN = TRUE,
                                             DATA_AUG = FALSE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                             PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

#PLOT 2
sim_vals10 = list(m1 = 1.1, m2 = 0.1, m3 = 10)
df_results_sseb = PLOT_SS_MCMC_GRID(epi_data_sseb, mcmc_sseb_output10, n_mcmc, sim_vals10,
                                    mcmc_specs = list(model_type = 'SSEB',
                                                      mod_par_names = c('alpha', 'beta', 'gamma'),
                                                      seed_count = seedX, burn_in_pc = 0.05, thinning_factor = 10),
                                    priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                                       c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                                    FLAGS_LIST = list(SSI = FALSE, BURN_IN = FALSE, THIN = TRUE,
                                                      DATA_AUG = FALSE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                                      PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
