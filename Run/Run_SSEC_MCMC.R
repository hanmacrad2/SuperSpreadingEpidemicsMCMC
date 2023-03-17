#RUN SSEB MCMC 
#FOLDERS
library(MASS)
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R'); source('R/EPI_FUNCTIONS.R')
source('R/SSEC_MCMC.R')

#OUTPUT FOLDER
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/ssec"
seedX = 4

#CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

#PARAMETERS 
n_mcmc = 100000 #50000 
#alphaX = 0.8; num_days = 50 #110 #Question: 50 days ok? Run 110 now. 1000 runs takes x mins

#1.SIMULATE DATA
seedX = seedX + 1
set.seed(seedX)
epi_data_nb = SIMULATE_EPI_SSEC()
plot.ts(epi_data_ssec)
#SAVE DATA
saveRDS(epi_data_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb_', seedX, '.rds' ))
epi_data_sseb = sim7

# *****************
#PARAM 1 FLAG_NEGBIN_PARAMATERISATION = list(param_prob = FALSE, param_mu = TRUE))
#*****************

#2.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_nb = MCMC_INFER_SSEC(epi_data_nb, n_mcmc = 100000)  #epi_data_ssec
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_nb$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_nb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssec_', seedX, '.rds' ))

#PARAMS
mcmc_nb_vec = mcmc_nb$ssec_params_matrix
k1_mcmc = mcmc_nb_vec[,1]; r01_mcmc = mcmc_nb_vec[,2]
#BURN IN
burn_in = 1000 #Equivalent to 10k as thinned 
k1_mcmc = k1_mcmc[burn_in:length(k1_mcmc)]

#CONFIDENCE INTERVALS
get_lower_ci(k1_mcmc); get_upper_ci(k1_mcmc)

get_lower_ci(r01_mcmc); get_upper_ci(r01_mcmc)

#PLOTS
par(mfrow = c(2,1))
plot.ts(k1[2000:length(k1)], main = 'k (dispersion parameter) - Negative Binomial model of daily infections. ', 
        ylab = 'k')
plot.ts(r01[2000:length(r01)], main = 'R0 - Negative Binomial model (Burn-in = 2000, mcmc reps = 100k)', ylab = 'R0')
plot.ts(mcmc_nb_vec)

#m1 mean
k1_mean = cumsum(k1)/seq_along(k1)
plot(seq_along(k1_mean), k1_mean, #ylim= m1_lim,
     xlab = 'Time', ylab =  'k1_mean',
     main = 'k1_mean', #bquote(bold(R[0] ~ "MCMC mean, Start:" ~ .(mcmc_specs$mod_start_points$m1))),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)


# *****************
#PARAM 2 
#*****************

#I. RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_nb2 = MCMC_INFER_SSEC(epi_data_nb, n_mcmc = 100000, FLAG_NEGBIN_PARAMATERISATION = list(param_prob = TRUE, param_mu = FALSE))  #epi_data_ssec
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_nb2$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_nb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssec_', seedX, '.rds' ))

#PARAMS
mcmc_nb_vec = mcmc_nb$ssec_params_matrix
k1_mcmc = mcmc_nb_vec[,1]; r01_mcmc = mcmc_nb_vec[,2]
#BURN IN
burn_in = 1000 #Equivalent to 10k as thinned 
k1_mcmc = k1_mcmc[burn_in:length(k1_mcmc)]

#CONFIDENCE INTERVALS
get_lower_ci(k1_mcmc); get_upper_ci(k1_mcmc)

get_lower_ci(r01_mcmc); get_upper_ci(r01_mcmc)

#SEED
seedX = seedX + 1