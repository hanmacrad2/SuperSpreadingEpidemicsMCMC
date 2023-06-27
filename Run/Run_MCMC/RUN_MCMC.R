#RUN MCMC
library(SuperSpreadingEpidemicsMCMC)
library(MASS)

#FOLDER STRUCTURE
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/"

OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/DATA_BASELINE_1/')

OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSEB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSNB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIR_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIB_DATA/DATA_SSIB_2/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIB_DATA/DATA_SSIB_2/')

OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIB_DATA/DATA_SSIB_4/')

OUTER_FOLDER = paste0(OUTER_FOLDER, 'MOCK_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/MISSING/')

#PARAMS
run = 1
n_repeats = 10; NMCMC = 50000
EPI_DATA = data_baseline

EPI_DATA = data_ssib2
EPI_DATA = data_ssib4

#***********************
# 1. RUN BASELINE MCMC
#**********************
run = '1_exp_1'
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE))

run = '1_ga'
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE))

#FLAGS_LIST = list(ADAPTIVE = TRUE, PRIOR_EXP = FALSE, PRIOR_GAMMA = TRUE, THIN = TRUE, BURN_IN = TRUE))


#***********************
# 3. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run,
                        n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = FALSE)) #data_sseb

#***********************
# 5. RUN SSIB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = 10, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = FALSE, SSIB = TRUE))
#run = 'run_2_ga_prior'

#***********************
# 2. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run,
                        n_repeats = n_repeats, n_mcmc = NMCMC,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIR = FALSE, SSIB = FALSE)) #data_ssnb

#**************************
# 4. RUN SSIR MCMC
#**************************
RUN_MCMC_MULTIPLE_TIMES(EPI_DATA, OUTER_FOLDER, run_number = run, n_repeats = 10, n_mcmc = 100000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIR = TRUE, SSIB = FALSE)) #data_ssir


#mcmc_ssib = MCMC_INFER_SSIB(data_ssib, n_mcmc = 1000)


#***********************
# 2. RUN SSNB MCMC
#**********************
run = 'gamma_prior_k'
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run,
                        n_repeats = 50, n_mcmc = 30000,
                        PRIORS_USED = list(EXP_K = FALSE, GAMMA_K = TRUE),
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIR = FALSE))