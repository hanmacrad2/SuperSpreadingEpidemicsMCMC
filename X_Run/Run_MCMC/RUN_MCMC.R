#RUN MCMC

#MCMC PARAMS
run = 1
n_repeats = 5; NMCMC = 30000
EPI_DATA = data_ssib

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
# 5. RUN SSI MCMC
#**************************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = TRUE, SSIB = FALSE)) 
