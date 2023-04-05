#NEW MCMC SCRIPT MULTI TIMES SCRIPT_2

#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSEB_DATA/')
create_folder(OUTER_FOLDER)

#***********************
# EPIDEMIC DATA (RUN AUTOMATICALLY)
#**********************
data_ssebI = SIMULATE_EPI_SSEB()
plot.ts(data_ssebI)
run_number = 1
saveRDS(data_ssebI, file = paste0(OUTER_FOLDER, 'epi_data_sseb_1.rds'))

#***********************
# 1. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_ssebI, OUTER_FOLDER, run_number = 1, n_repeats = 500, n_mcmc = 500000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_ssebI, OUTER_FOLDER, run_number = 1, n_repeats = 500, n_mcmc = 500000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_ssebI, OUTER_FOLDER, run_number = 1, n_repeats = 500, n_mcmc = 500000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))
