#IMPLEMENT MODEL EVIDENCE VIA IMPORTANCE SAMPLING 
library(SuperSpreadingEpidemicsMCMC)

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#***********************
# 1. BASE DATA (RUN AUTOMATICALLY)
#**********************

#1. Data
BASE_DATA_LOC = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(BASE_DATA_LOC, 'epi_data_base_1.rds'))

#FOLDER SAVE (SSIB)
OUTPUT_FOLDER = paste0(BASE_DATA_LOC, 'SSIB/')

#SSIB(SSIB = TRUE)
RUN_MCMC_MODEL_EV_IMP_SAMP(data_baseI, OUTPUT_FOLDER)

#***************************
# 2. MANUAL ITERATION
#***************************

#ITERATION
run = 1
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)
create_folder(CURRENT_OUTPUT_FOLDER)

#***********
#* 1. BASE
#***********

ests_base = c()
for (i in 1:100){
  print(paste0('i = ', i))
  
  #1. MCMC samples
  mcmc_base = MCMC_INFER_BASELINE(data_baseI)
  saveRDS(mcmc_base, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base_', i ))
  #2. Model 
  phat_base = GET_IMP_SAMP_MODEL_EV_BASE(mcmc_base$r0_vec, data_baseI) #NA
  ests_base[i] = phat_base
  print(ests_base)
}


#***********
#* 2. SSEB
#* **********
#FOLDER SAVE
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/SSEB"

#ITERATION
run = 1
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)
create_folder(CURRENT_OUTPUT_FOLDER)

ests_sseb = c()

for (i in 1:100){
  print(paste0('i = ', i))
  
  #1. MCMC samples
  n_mcmc = 30000
  mcmc_sseb = MCMC_INFER_SSEB(data_baseI, n_mcmc)
  saveRDS(mcmc_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', i ))
  #2. Model 
  mcmc_samples =  matrix(c(mcmc_sseb$alpha_vec, mcmc_sseb$beta_vec, mcmc_sseb$gamma_vec), ncol = 3)
  
  phat_sseb = GET_IMP_SAMP_MODEL_EV_SSEB(mcmc_samples, data_baseI) #NA
  ests_sseb[i] = phat_sseb
  print(ests_sseb)
}


#RUNNING ABOVE

# *******************************************************
# PLOT
#*******************************************************


#*********************************
#1. MODEL EVIDENCE BASELINE
#r0_samples = c(1.2, 1.3, 1.31, 1.21, 1.25, 1.22, 1.23, 1.24)
#ests_base = c()

for (i in 1:100){
  print(paste0('i = ', i))
  phat_base = GET_IMP_SAMP_MODEL_EV_BASE(r0_samples, data_baseI) #NA
  ests_base[i] = phat_base
  print(ests_base)
}

#PLOTS
par(mfrow = c(2,1))
plot(seq_along(ests_base), ests_base,
     ylim = c(min(ests_base)-50, max(ests_base)+50), lwd = 1, pch = 19,
     main = 'Model evidence (Baseline) via Importance sampling. Base data',
     xlab = 'rep', ylab = 'Model evidence estimate')

#PLOT 2 (Close up)
plot(seq_along(ests_base), ests_base,
     ylim = c(min(ests_base)-0.5, max(ests_base)+0.5), lwd = 1, pch = 19,
     main = 'Model evidence (Baseline) via Importance sampling. Base data. (Close up)',
     xlab = 'rep', ylab = 'Model evidence estimate')


#************************************************
#2. RUN SSEB
ests = c()
ests_sseb = ests

for (i in 1:100){
  print(paste0('i = ', i))
  phat1 = GET_IMP_SAMP_MODEL_EV_SSEB(mcmc_samples, data_baseI) #NA
  ests[i] = phat1
  print(ests)
}

#PLOTS
par(mfrow = c(2,1))
plot(seq_along(ests), ests,
     ylim = c(min(ests)-50, max(ests)+50), lwd = 1, pch = 19,
     main = 'Model evidence (SSEB) via Importance sampling. Base data',
     xlab = 'rep', ylab = 'Model evidence estimate')

#PLOT 2 (Close up)
plot(seq_along(ests), ests,
     ylim = c(min(ests)-0.5, max(ests)+0.5), lwd = 1, pch = 19,
     main = 'Model evidence (SSEB) via Importance sampling. Base data. (Close up)',
     xlab = 'rep', ylab = 'Model evidence estimate')


#***************************
#* APPLY FUNCTION WITH TOY EXAMPLE
#***************************

alpha = c(0.8, 0.81, 0.805, 0.9, 0.10)
beta = c(0.02, 0.025, 0.03, 0.5, 0.6)
gamma = c(10, 10.2, 10.1, 11, 12)
mcmc_samples = matrix(c(alpha, beta, gamma), ncol = 3)

rlnorm.rplus(100,log(mean(mcmc_samples)),cov(mcmc_samples))
dlnorm.rplus(x,log(mean(mcmc_samples)), cov(mcmc_samples))
imp_samp_comps1 = GET_PROPOSAL_MULTI_DIM(mcmc_samples, data_baseI) 

#ESTIMATE
phat1 = GET_IMP_SAMP_MODEL_EV_SSEB(mcmc_samples, data_baseI) #NA
phat1

phat_sseb = GET_IMP_SAMP_MODEL_EV_SSEB(mcmc_samples, data_sseb1) #NA
phat_sseb

#RUN MULTIPLE TIMES

#**********
#* OTHER

#PLOT
par(mfrow = c(1,1))

plot.ts(ests,
        xlab = 'rep', ylab = 'Estimate',
        main = 'Model Evidence estimate, repeat (same dataset). lognormal propasls = variable. Use same?')

plot(seq_along(ests), ests,
     xlab = 'rep', ylab = 'Estimate',
     main = 'Model Evidence estimate, repeat (same dataset). lognormal propasls = variable. Use same?')


#SSIB RESULTS
ests_ssib = c(-17.97085, -17.02987, -17.66838, -17.88493, -16.89313, -17.44772, -18.76960, -17.65383,
              -18.59628, -17.03140, -16.74319, -17.76785, -17.46280, -17.16593, -18.31888, -17.84318,
              -16.95889, -17.51474, -17.38049, -17.59530, -18.18540, -17.35371, -17.78427, -18.00840,
              -17.36284, -17.00414, -17.09351, -17.34403, -17.41790, -17.21114, -17.16511, -17.66789,
              -17.38434, -17.56134, -17.05074, -17.51592, -17.91808, -17.35481, -17.10992, -17.41607,
              -17.33731, -17.81293, -18.14597, -18.11503, -19.79129, -19.54160, -18.28430, -17.42519,
              -17.45460, -17.04829, -17.64163, -17.90642, -17.70591, -16.89735, -17.10859,
              -17.22414, -17.77228, -16.78630, -16.95084, -18.34672, -17.14331, -17.62390, -16.92714,
              -17.66143, -16.55392, -18.07520, -16.94167, -17.26395, -18.13378, -17.43426, -17.47782,
              -16.76090, -17.31100, -17.39757, -17.44080, -17.05849, -17.51593, -16.68125, -17.28342,
              -16.94697, -17.98902, -17.25548, -18.30939, -17.34057, -16.47667, -17.31498, -17.95910,
              -17.45501, -18.58359, -17.42545, -17.14935, -17.72494, -17.63773, -17.42719, -16.89223,
              -17.77827, -17.11730, -18.23225, -16.99389, -16.99354)
