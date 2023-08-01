#RUN MCMC

#MCMC PARAMS
run = 1
n_repeats = 5; NMCMC = 30000
EPI_DATA = data_ssnb2

#***********************
# 1. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE))

#***********************
# 2. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run,
                        n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIR = FALSE, SSIB = FALSE)) #data_ssnb

#***********************
# 3. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run,
                        n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE)) #data_sseb

#***********************
# 4. RUN SSIB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = TRUE))

#**************************
# 5. RUN SSIR MCMC
#**************************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = TRUE, SSIB = FALSE)) 

#******************************************************************************
#*
#* OTHER PRIORS
#* 
#*****************************************************************************

#run = '1_exp_1'
#run = '1_ga'; #run = 'run_2_ga_prior'
# RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
#                         FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
#                                             SSIR = FALSE, SSIB = FALSE))

#***********************
# 2. RUN SSNB MCMC
#**********************
run = 'gamma_prior_k'
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run,
                        n_repeats = 50, n_mcmc = 30000,
                        PRIORS_USED = list(EXP_K = FALSE, GAMMA_K = TRUE),
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIR = FALSE))