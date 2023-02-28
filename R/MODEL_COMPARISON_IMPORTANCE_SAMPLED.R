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
#1. GET PROPSAL
GET_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data, #priors = 
                                   n_samples = 100) {               #1.OTHER: GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #DEFENSE MIXTURE 
  means = colMeans(mcmc_samples)
  theta_samples = rlnorm.rplus(n_samples, log(colMeans(mcmc_samples)), cov(mcmc_samples))
  #print(paste0('theta_samples = ', theta_samples))
  proposal = dlnorm.rplus(theta_samples, log(colMeans(mcmc_samples)), cov(mcmc_samples))
  #print(paste0('proposal = ', proposal))
  
  imp_samp_comps = list(theta_samples = theta_samples, proposal = proposal)
  
  #t dist (incompatible as not strictly positive)
  #theta_samples = rmvt(n_samples, sigma = cov(mcmc_samples), df = 3) + means #Samples (theta)
  #proposal = dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = FALSE)
  
  #PRIORS
  #priors = dexp(mcmc_samples[,1]) + dexp(mcmc_samples[,2]) + dexp((mcmc_samples[,3] - 1)) #LOGS
  #MIXTURE
  #q_defense_mixture = 0.95*proposal + 0.05*priors
  
  return(imp_samp_comps)
}

#**********
#* PHAT SSEB MODEL
GET_IMP_SAMP_MODEL_EV_SSEB <-function(mcmc_samples, epidemic_data, n_samples = 100) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  imp_samp_comps = GET_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data)
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


# if (is.na(estimate)){
#   print(paste0('loglikesseb', LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
#                                             theta_samples[i, 3])))
#   print(paste0(' log(priors[i])',  log(priors[i])))
#   print(paste0(' log(q_defense_mixture[i])', log(q_defense_mixture[i])))
# }