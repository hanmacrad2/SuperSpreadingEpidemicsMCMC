#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASE_DATA/')

#***********************
# EPIDEMIC DATA (RUN AUTOMATICALLY)
#**********************

# BASE DATA
data_baseI = readRDS(file = paste0(OUTER_FOLDER, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
run_number = 1

#2. SSEB DATA 
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSEB_DATA/')
model_type = 'sseb'; run = 1
data_sseb = readRDS(file = paste0(OUTER_FOLDER, 'epi_data_', model_type, '_', run, '.rds'))
plot.ts(data_sseb)

#***********************
# 1. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_baseI, OUTER_FOLDER) 

#***********************
# 2. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_baseI, OUTER_FOLDER, run_number = 2, n_repeats = 500, n_mcmc = 500000,
                                    FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                        SSIB = FALSE, SSIC = FALSE))

#RUN_MULTIPLE_MCMC_SSEB(data_baseI, OUTER_FOLDER, n_reps = 100, n_mcmc = 100000) 

#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_baseI, OUTER_FOLDER, run_number = 2, n_repeats = 500, n_mcmc = 500000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))

#RUN_MULTIPLE_MCMC_SSNB(data_baseI, OUTER_FOLDER, n_reps = 100, n_mcmc = 100000) 

#***********************
# 4. RUN SSIC MCMC
#**********************
RUN_MULTIPLE_MCMC_SSIC(data_baseI, OUTER_FOLDER, n_reps = 100, n_mcmc = 100000) 

#**********************
#Inspect MCMC results
#**********************
model_type = 'SSEB'; i = 100; run_number = 2
RESULTS_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')

#MCMC RESULTS
mcmc50 = readRDS(file = paste0(RESULTS_FOLDER, 'mcmc_', model_type , '_', i, '.rds'))
PLOT_SSB_MCMC_REAL_DATA(data_sseb, mcmc50, 50000)

PLOT_SSB_MCMC_GRID(data_sseb, mcmc50, 50000)

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

