#NEW MCMC SCRIPT MULTI TIMES SCRIPT_3

#RUN MULTIPLE MCMC ITERATIONS
library(SuperSpreadingEpidemicsMCMC)

#FOLDER RESULTS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'NZ_DATA_WAIT_21/')
create_folder(OUTER_FOLDER)
run_number = 2

#***********************
# EPIDEMIC DATA -- NEW ZEALAND
#**********************
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
data_file_wait_21 = read.csv(paste0(DATA_FOLDER, 'data_waitemata_aug_21.csv'))
data_wait_08_21 = data_file_wait_21$Cases
plot.ts(data_wait_08_21, ylab = 'Infection count', main = 'Waitemata NZ, August 2021')
saveRDS(data_wait_08_21, file = paste0(OUTER_FOLDER, 'data_wait_08_21.rds'))

#***********************
# 2. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 20000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 20000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 20000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))


#INSPECT
mcmc_output = MCMC_INFER_BASELINE(data_wait_08_21, n_mcmc = 1000)

mcmc_output = MCMC_INFER_SSEB(data_wait_08_21, n_mcmc = 1000)   


for (t in 2:num_days) {
  
  #ETA (t-1)
  eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = R0X/k) #Draw eta from previous time step
  #INFECTIVITY
  infectivity = rev(prob_infect[1:t]) 
  #POISSON; OFFSPRINT DISTRIBUTION
  total_rate = sum(eta_vec*infectivity) #DOT PRODUCT
  x[t] = rpois(1, total_rate)
  
}
