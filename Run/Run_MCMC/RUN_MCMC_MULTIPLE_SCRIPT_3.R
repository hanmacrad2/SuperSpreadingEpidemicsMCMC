#NEW MCMC SCRIPT MULTI TIMES SCRIPT_3

#RUN MULTIPLE MCMC ITERATIONS
library(SuperSpreadingEpidemicsMCMC)

DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'NZ_DATA_WAIT_21/')
create_folder(OUTER_FOLDER)
run_number = 1

#***********************
# EPIDEMIC DATA -- NEW ZEALAND
#**********************
data_file_wait_21 = read.csv(paste0(DATA_FOLDER, 'data_waitemata_aug_21.csv'))
data_wait_08_21 = data_file_wait_21$Cases
plot.ts(data_wait_08_21)
saveRDS(data_wait_08_21, file = paste0(OUTER_FOLDER, 'data_wait_08_21.rds'))

#***********************
# 2. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = 1, n_repeats = 100, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = 1, n_repeats = 100, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = 1, n_repeats = 100, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))
