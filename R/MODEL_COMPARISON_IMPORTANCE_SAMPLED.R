#Model Comparison - Importance Sampled
#'Model evidence estimator via importance sampling'

#LIBRARIES
library(SuperSpreadingEpidemicsMCMC)
library(compositions)
#library(mvtnorm)

#**********
#*LOG LIKELIHOOD
LOG_LIKE_SSEB <- function(x, lambda_vec, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x); logl = 0
  
  for (t in 2:num_days) {
    
    #Time
    #print(paste0('time =', t))
    #lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum_xt = 0
    term1 = exp(-alphaX*lambda_vec[t]); term2 = alphaX*lambda_vec[t]
    #print(paste0('term1 =', term1)); print(paste0('term2 =', term2))
    
    for (nt in 0:x[t]){ #Sum for all possible values of nt, xt-nt
      
      #Log likelihood
      st = x[t] - nt
      prob_st =  PROBABILITY_ST(st, lambda_vec[t], alphaX, betaX, gammaX)
      if (is.na(prob_st)){
        print(paste0('prob_st  =', prob_st))
        print(paste0('nt  =', nt))
        print(paste0('alphaX  =', alphaX))
        print(paste0('betaX  =', betaX))
        print(paste0('gammaX  =', gammaX))
      }
      
      inner_sum_xt = inner_sum_xt + 
        term1*(term2)^nt*(1/factorial(nt))*prob_st
    } 
    #print(paste0('log(inner_sum_xt)  =', log(inner_sum_xt) ))
    logl = logl + log(inner_sum_xt) 
  }
  
  return(logl)
}

#*****************
GET_PROPOSAL_UNI_VAR <- function(mcmc_samples, epidemic_data, #priors = 
                                   n_samples) {               #1.OTHER: GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #PROPOSAL 
  means = mean(mcmc_samples)
  theta_samples = rlnorm(n_samples, log(mean(mcmc_samples)), sd(mcmc_samples))
  proposal = dlnorm(theta_samples, log(mean(mcmc_samples)), sd(mcmc_samples))
  imp_samp_comps = list(theta_samples = theta_samples, proposal = proposal)
  
  return(imp_samp_comps)
}

#***************************************
#1B. GET IMPORTANCE SAMPLING PROPOSAL
#***************************************
GET_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data, #priors = 
                                   n_samples) {               #1.OTHER: GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #PROPOSAL 
  means = colMeans(mcmc_samples)
  theta_samples = rlnorm.rplus(n_samples, log(colMeans(mcmc_samples)), cov(mcmc_samples))
  proposal = dlnorm.rplus(theta_samples, log(colMeans(mcmc_samples)), cov(mcmc_samples))
  imp_samp_comps = list(theta_samples = theta_samples, proposal = proposal)
  
  #t dist (incompatible as not strictly positive)
  #theta_samples = rmvt(n_samples, sigma = cov(mcmc_samples), df = 3) + means #Samples (theta)
  #proposal = dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = FALSE)
  
  return(imp_samp_comps)
}


#***************************************************
#* GET MODEL EVIDENCE ESTIMATES VIA IMPORTANCE SAMPLING
#*****************************************************

#**********************************************
#* PHAT SSEB MODEL
GET_IMP_SAMP_MODEL_EV_BASE <-function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0
  imp_samp_comps = GET_PROPOSAL_UNI_VAR(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  proposal = imp_samp_comps$proposal
  
  #PRIORS 
  priors = dexp(theta_samples)  #dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    
    estimate = LOG_LIKE_BASELINE(epidemic_data, theta_samples[i]) +
      log(priors[i]) - log(q_defense_mixture[i])
    
    if(is.infinite(estimate)){
      print('yes infinite')
      print(paste0('LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])', LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])))
      print(paste0('log(priors[i])', log(priors[i])))
      print(paste0('theta_samples[i]', theta_samples[i]))
      print(paste0('log(q_defense_mixture[i])', log(q_defense_mixture[i])))
    }
    
    #print(estimate)
    sum_estimate = sum_estimate + estimate
  }
  
  p_hat_est = sum_estimate/n_samples
  
  return(p_hat_est)
}

#**********************************************
#* PHAT SSEB MODEL
GET_IMP_SAMP_MODEL_EV_SSEB <-function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  imp_samp_comps = GET_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  proposal = imp_samp_comps$proposal
  
  #PRIORS 
  priors = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    
    estimate = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
                             theta_samples[i, 3]) +
      log(priors[i]) - log(q_defense_mixture[i])

    #print(estimate)
    sum_estimate = sum_estimate + estimate
  }
  
  p_hat_est = sum_estimate/n_samples
  
  return(p_hat_est)
}

