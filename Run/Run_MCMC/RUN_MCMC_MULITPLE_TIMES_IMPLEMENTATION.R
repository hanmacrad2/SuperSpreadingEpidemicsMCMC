#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#***********************
# 1. EPIDEMIC DATA (RUN AUTOMATICALLY)
#**********************
LOC_BASE_DATA = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(LOC_BASE_DATA, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
run_number = 1

#***********************
# 2. RUN SSNB MCMC
#**********************
RUN_MULTIPLE_MCMC_SSNB(data_baseI, OUTER_FOLDER, n_reps = 100, n_mcmc = 100000) 

#Inspect MCMC results
model_type = 'SSNB'; i = 100
RESULTS_FOLDER = paste0(OUTER_FOLDER, '/', model_type, '/run_', run_number, '/')

#MCMC RESULTS
i = 90
mcmc90 = readRDS(file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
#PLOT
PLOT_SSNB_MCMC_GRID(data_baseI, mcmc90, n_mcmc = 100000)

ssnb1 = mcmc1$ssec_params_matrix
plot.ts(ssec1)

#BURN IN
k1 = ssec1[,1]; r01 = ssec1[,2]
par(mfrow = c(2,1))
plot.ts(k1[2000:length(k1)], main = 'k (dispersion parameter) - Negative Binomial model of daily infections. ', 
        ylab = 'k')
plot.ts(r01[2000:length(r01)], main = 'R0 - Negative Binomial model (Burn-in = 2000, mcmc reps = 100k)', ylab = 'R0')

i = 5 
mcmc5 = readRDS(file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
ssec5 = mcmc5$ssec_params_matrix
plot.ts(ssec5)

mcmc100 = readRDS(file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
ssec100 = mcmc100$ssec_params_matrix
plot.ts(ssec100)
