#NEW MCMC SCRIPT MULTI TIMES SCRIPT_2

#RUN MULTIPLE MCMC ITERATIONS
#library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#MODEL TYPE
model_type = 'sseb';
OUTER_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '_DATA/')
create_folder(OUTER_FOLDER)
run_number = 2

#***********************
# EPIDEMIC DATA (RUN AUTOMATICALLY)
#**********************
data_ssebI = SIMULATE_EPI_SSEB()
plot.ts(data_ssebI)
data_num = 1
saveRDS(data_sseb, file = paste0(OUTER_FOLDER, 'epi_data_', model_type, '_', data_num, '.rds'))

#LOAD DATA
data_sseb = readRDS(file = paste0(OUTER_FOLDER, 'epi_data_', model_type, '_', data_num, '.rds'))
plot.ts(data_sseb)

#***********************
# 1. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_sseb, OUTER_FOLDER, run_number = run_number, n_repeats = 100, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_sseb, OUTER_FOLDER, run_number = run_number, n_repeats = 100, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))


#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_sseb, OUTER_FOLDER, run_number = run_number, n_repeats = 100, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))
