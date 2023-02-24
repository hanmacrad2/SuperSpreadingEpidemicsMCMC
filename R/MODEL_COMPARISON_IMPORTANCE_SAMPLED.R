#Model Comparison - Importance Sampled

#'Model evidence estimator via importance sampling'

#LIBRARIES
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)

#*****************
#1. GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
GET_MULTI_DIM_PROPOSAL <- function(mcmc_samples, epidemic_data, n_samples = 100) {
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #DEFENSE MIXTURE 
  means = colMeans(mcmc_samples)
  theta_samples = rmvt(n_samples, sigma = var(mcmc_samples), df = 3) + means #Samples (theta)
  proposal = dmvt(theta_samples - means, sigma = var(mcmc_output), df = 3, log = FALSE)
  
  #PRIORS
  priors = dexp(mcmc_samples[,1]) + dexp(mcmc_samples[,2]) + dexp((mcmc_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  return(proposal)
}

#1. GET PROPSAL
GET_MULTI_DIM_PROPOSAL <- function(mcmc_samples, epidemic_data, n_samples = 100) {
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #DEFENSE MIXTURE 
  means = colMeans(mcmc_samples)
  theta_samples = rmvt(n_samples, sigma = cov(mcmc_samples), df = 3) + means #Samples (theta)
  proposal = dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = FALSE)
  
  #PRIORS
  priors = dexp(mcmc_samples[,1]) + dexp(mcmc_samples[,2]) + dexp((mcmc_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  return(proposal)
}

#**********
#* PHAT SSEB MODEL
GET_MODEL_EVIDENCE_EST_BASELINE <-function(mcmc_samples, proposal, epidemic_data, n_samples = 100) {
  
  'Estimate of model evidence using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #PRIORS exp(1)
  priors = dexp(mcmc_samples[,1])
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    sum_estimate = sum_estimate +
      LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],  theta_samples[i, 3]) +
      log(priors[i]) - log(q_defense_mixture[i])
  }
  
  p_hat_est = sum_estimate/n_samples
}

#**********
#* PHAT SSEB MODEL
GET_MODEL_EVIDENCE_EST_SSEB <-function(mcmc_samples, proposal, epidemic_data, n_samples = 100) {
  
  'Estimate of model evidence using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #PRIORS
  priors = dexp(mcmc_samples[,1]) + dexp(mcmc_samples[,2]) + dexp((mcmc_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    sum_estimate = sum_estimate +
      LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],  theta_samples[i, 3]) +
      log(priors[i]) - log(q_defense_mixture[i])
  }
  
  p_hat_est = sum_estimate/n_samples
}

#INSPECT LIKELIHOOD
LOG_LIKE_SSEB <- function(x, lambda_vec, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x); logl = 0
  
  for (t in 2:num_days) {
    
    #Time
    print(paste0('time =', t))
    #lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum_xt = 0
    term1 = exp(-alphaX*lambda_vec[t]); term2 = alphaX*lambda_vec[t]
    #print(paste0('term1 =', term1)); print(paste0('term2 =', term2))
    
    for (nt in 0:x[t]){ #Sum for all possible values of nt, xt-nt
      
      #Log likelihood
      st = x[t] - nt
      prob_st =  PROBABILITY_ST(st, lambda_vec[t], alphaX, betaX, gammaX)
      print(paste0('prob_st  =', prob_st))
      
      inner_sum_xt = inner_sum_xt + 
        term1*(term2)^nt*(1/factorial(nt))*prob_st
    } 
    print(paste0('log(inner_sum_xt)  =', log(inner_sum_xt) ))
    logl = logl + log(inner_sum_xt) 
  }
  
  return(logl)
}

#************************************
#* 2. APPLY WITH MCMC OUTPUT BASELINE
#************************************
#Do multiple runs over the weekend 

#***************************
#* APPLY FUNCTION WITH TOY EXAMPLE
#***************************

alpha = c(0.8, 0.81, 0.805, 0.9, 0.10)
beta = c(0.02, 0.025, 0.03, 0.5, 0.6)
gamma = c(10, 10.2, 10.1, 11, 12)
mcmc_samples = matrix(c(alpha, beta, gamma), ncol = 3)
lambda_vec = get_lambda(data_baseI)

phat1 = ESTIMATE_MODEL_EVIDENCE(mcmc_samples, data_baseI) #NA
phat1





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

#TOTAL FUNCTION OF ESTIMATE OF MODEL EVIDENCE
ESTIMATE_MODEL_EVIDENCE <-function(mcmc_samples, epidemic_data, n_samples = 100) {
  
  'Estimate of model evidence using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #DEFENSE MIXTURE 
  means = colMeans(mcmc_samples)
  theta_samples = rmvt(n_samples, sigma = cov(mcmc_samples), df = 3) + means #Samples (theta)
  proposal = dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = FALSE)
  
  #PRIORS
  priors = dexp(mcmc_samples[,1]) + dexp(mcmc_samples[,2]) + dexp((mcmc_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    sum_estimate = sum_estimate +
      LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],  theta_samples[i, 3]) +
      log(priors[i]) - log(q_defense_mixture[i])
  }
  
  p_hat_est = sum_estimate/n_samples
}
