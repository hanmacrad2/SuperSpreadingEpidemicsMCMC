#Model Comparison - Importance Sampled

#LIBRARIES
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)

#*****************
ESTIMATE_MODEL_EVIDENCE <-function(mcmc_samples, epidemic_data, n_samples = 100) {
  
  'Estimate of model evidence using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #DEFENSE MIXTURE 
  means = colMeans(mcmc_samples)
  proposal_samples = rmvt(n_samples, sigma = cov(mcmc_samples), df = 3) + means #Samples (theta)
  log_proposal = dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = TRUE) #FALSE #Should log == false? Maybe doesnt' matter which step. think about it.
  
  #PRIORS
  log_priors = dexp(mcmc_samples[,1], log = TRUE) + dexp(mcmc_samples[,2], log = TRUE) + dexp((mcmc_samples[,3] - 1), log = TRUE) #LOGS
  
  #MIXTURE
  log_q_defense_mixture = 0.95*log_proposal + 0.05*log_priors
  
  #LOOP
  for(i in 1:n_samples){
    sum_estimate = sum_estimate +
LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],  theta_samples[i, 3]) +
      log_priors[i] - log_q_defense_mixture[i]
  }
  
  p_hat_est = sum_estimate/n_samples
}

#******************************
#* 1. APPLY WITH MCMC OUTPUT BASELINE
#******************************


#******************
#* 2. APPLY WITH MCMC OUTPUT BASELINE
#******************
#Do multiple runs over the weekend 


#******************
#* APPLY FUNCTION WITH TOY EXAMPLE
#******************

#Trial
alpha = c(0.8, 0.81, 0.805, 0.9, 0.10)
beta = c(0.02, 0.025, 0.03, 0.5, 0.6)
gamma = c(10, 10.2, 10.1, 11, 12)
mcmc_samples = matrix(c(alpha, beta, gamma), ncol = 3)
lambda_vec = get_lambda(data_baseI)

phat1 = ESTIMATE_MODEL_EVIDENCE(mcmc_samples, data_baseI)

#******************
#* TOY EXAMPLE
#******************

#Trial
alpha = c(0.8, 0.81, 0.805, 0.9, 0.10)
beta = c(0.02, 0.025, 0.03, 0.5, 0.6)
gamma = c(10, 10.2, 10.1, 11, 12)
mcmc_output = matrix(c(alpha, beta, gamma), ncol = 3)
#Mean of mcmc
means = colMeans(mcmc_output)

#DEFENSE MIXTURE #Should log == false?
theta_samples = rmvt(1000, sigma = cov(mcmc_output), df = 3) + means #Samples (theta)
q_defense = 0.95*dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = FALSE)
#Priors
priors = 0.5*dexp(alpha) + 0.5*dexp(beta) + 0.5*dexp(gamma - 1)
#Q mixture defense
q_defense = q_defense + priors

#P Hat estimator

get_p_hat_estimator <- function(q_defense_mixture){
  
  #Loop over sample
  sum_estimate = 0
  
  for(i in 1:N){
    q_defense_mixture[i]
  }
}

get_p_hat_estimator <- function(q_defense_mixture){
  
  #Loop over sample
  sum_estimate = 0
  N = legnth(q_defense_mixture)
  
  for(i in 1:N){
    log(q_defense_mixture[i]) + LOG_LIKE_SSEB(q_defense_mixture[i]) - q_defense_mixture[i]
  }
}



#Estimator for one model
#Do for each model :) 


#LOAD MCMC
