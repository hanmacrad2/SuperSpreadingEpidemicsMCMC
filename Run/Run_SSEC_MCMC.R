#RUN SSEB MCMC 
#FOLDERS
library(MASS)
library(SuperSpreadingEpidemicsMCMC)

#OUTPUT FOLDER
seedX = 1; n_mcmc = 100000
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/SSNB"
FOLDERX = paste0(OUTPUT_FOLDER,  '/run_', seedX)
create_folder(FOLDERX)

#1.SIMULATE DATA
params_nb = c(0.16, 1.2)
seedX = seedX + 1; set.seed(seedX)
epi_data_nb = SIMULATE_EPI_SSNB()
plot.ts(epi_data_nb)
#Likelihood
lambda_vec1 = get_lambda(epi_data_nb)
likenb1 = LOG_LIKE_SSNB(epi_data_nb, lambda_vec1, params_nb) 
#SAVE DATA
saveRDS(epi_data_nb, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_ssnb_', seedX, '.rds' ))

#*********************************************
#PARAMERISATION 1: NEGATIVE BINOMIAL MODEL (k, mean)
#*********************************************
#i.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_nb = MCMC_INFER_SSNB(epi_data_nb, n_mcmc = 100000,
                          FLAG_NEGBIN_PARAMATERISATION = list(param_prob = FALSE, param_mu = TRUE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_nb$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_nb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_SSNB_', seedX, '.rds' ))

#PARAMS
mcmc_nb_vec = mcmc_nb$ssnb_params_matrix
k1_mcmc = mcmc_nb_vec[,1]; r01_mcmc = mcmc_nb_vec[,2]
#BURN IN
burn_in = 1000 #Equivalent to 10k as thinned 
k1_mcmc = k1_mcmc[burn_in:length(k1_mcmc)]
r01_mcmc = r01_mcmc[burn_in:length(r01_mcmc)]

#PLOT
PLOT_SSNB_MCMC_GRID(epi_data_nb, mcmc_nb, n_mcmc, likenb1)
                                
#CONFIDENCE INTERVALS
get_lower_ci(k1_mcmc); get_upper_ci(k1_mcmc)
get_lower_ci(r01_mcmc); get_upper_ci(r01_mcmc)

#PLOTS
par(mfrow = c(2,1))
plot.ts(k1, main = 'k (dispersion parameter) - Negative Binomial model of daily infections. ', 
        ylab = 'k')
plot.ts(k1, main = 'R0 - Negative Binomial model (Burn-in = 2000, mcmc reps = 100k)', ylab = 'R0')
plot.ts(mcmc_nb_vec)

#PLOT CUM MEAN
PLOT_CUM_MEAN(k1_mcmc)

#*****************
#PARAMETERISATION 2 
#*****************

#I.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_nb2 = MCMC_INFER_SSNB(epi_data_nb, n_mcmc = 100000, FLAG_NEGBIN_PARAMATERISATION = list(param_prob = TRUE, param_mu = FALSE))  #epi_data_SSNB
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_nb2$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_nb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_SSNB_', seedX, '.rds' ))

#PARAMS
mcmc_nb_vec2 = mcmc_nb2$SSNB_params_matrix
k2_mcmc = mcmc_nb_vec2[,1]; r02_mcmc = mcmc_nb_vec2[,2]
#BURN IN
burn_in = 1000 #Equivalent to 10k as thinned 
k2_mcmc = k2_mcmc[burn_in:length(k2_mcmc)]; r02_mcmc = r02_mcmc[burn_in:length(r02_mcmc)]
#PLOT
plot.ts(k2_mcmc)
plot.ts(r02_mcmc)

#CONFIDENCE INTERVALS
get_lower_ci(k2_mcmc); get_upper_ci(k2_mcmc)
get_lower_ci(r02_mcmc); get_upper_ci(r02_mcmc)

#************************
#PLOT TOTAL RESULTS
par(mfrow = c(2,2))
plot.ts(k1_mcmc,
        main = 'k: Neg Bin ~ I(t). Params: (k, mean) ', 
        ylab = 'k')
plot.ts(k2_mcmc,
        main = 'k: Neg Bin ~ I(t). Params: (k, prob) ', 
        ylab = 'k')
plot.ts(r01_mcmc,
        main = 'R0: Neg Bin ~ I(t). Params: (k, mean) ', 
        ylab = 'R0')
plot.ts(r02_mcmc,
        main = 'R0: Neg Bin ~ I(t). Params: (k, prob) ', 
        ylab = 'R0')

#'k (dispersion parameter) - Negative Binomial model of daily infections. '
#'R0 - Negative Binomial model (Burn-in = 2000, mcmc reps = 100k)', ylab = 'R0'

#SEED
seedX = seedX + 1