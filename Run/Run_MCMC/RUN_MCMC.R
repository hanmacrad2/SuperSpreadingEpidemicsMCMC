#RUN MCMC
run = 1
library(MASS)
#FOLDER STRUCTURE
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/"

OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSEB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSNB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIR_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIB_DATA/')

#***********************
# 1. RUN BASELINE MCMC
#**********************
run = 1
RUN_MCMC_MULTIPLE_TIMES(data_baseline, OUTER_FOLDER, run_number = run, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE))


#***********************
# 2. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_ssnb, OUTER_FOLDER, run_number = run,
                        n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIR = FALSE, SSIB = FALSE))

#***********************
# 3. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE))

RUN_MCMC_MULTIPLE_TIMES(data_sseb, OUTER_FOLDER, run_number = run, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE))

#***********************
# 5. RUN SSIB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_ssib, OUTER_FOLDER, run_number = run, n_repeats = 10, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = TRUE))

#**************************
# 4. RUN SSIR MCMC
#**************************
RUN_MCMC_MULTIPLE_TIMES(data_ssir2$epidemic_data, OUTER_FOLDER, run_number = run, n_repeats = 10, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = TRUE, SSIB = FALSE))


#***********************
# 2. RUN SSNB MCMC
#**********************
run = 'gamma_prior_k'
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run,
                        n_repeats = 50, n_mcmc = 30000,
                        PRIORS_USED = list(EXP_K = FALSE, GAMMA_K = TRUE),
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIR = FALSE))