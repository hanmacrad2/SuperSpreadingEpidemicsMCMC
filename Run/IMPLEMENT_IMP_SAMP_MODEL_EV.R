#IMPLEMENT MODEL EVIDENCE VIA IMPORTANCE SAMPLING 
library(SuperSpreadingEpidemicsMCMC)

#FOLDER SAVE
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/BASE"

#ITERATION
run = 1
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)
create_folder(CURRENT_OUTPUT_FOLDER)

#***********
#* 1. BASE
#* **********

#TO DO
#1. GET MCMC SAMPLES AND CALCULATE ESTIMATES

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

#TO DO
#1. GET MCMC SAMPLES AND CALCULATE ESTIMATES

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


#CHEC
if(is.infinite(phat_base)){
  print('yes')
}