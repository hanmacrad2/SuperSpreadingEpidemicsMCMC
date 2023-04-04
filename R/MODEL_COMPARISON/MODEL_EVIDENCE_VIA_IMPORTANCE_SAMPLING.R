#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
library(compositions)

#*************************************
#* LOG SUM EXP
#* ***********************************
LOG_SUM_EXP <- function(vectorX){
  
  #REMOVE NA VALUES
  vectorX = na.omit(vectorX)
  
  max_val = max(vectorX)
  
  out = max_val + log(sum(exp(vectorX - max_val)))
  
  return(out)
}

#***************************************
#*
#1. GET IMPORTANCE SAMPLING PROPOSAL (LOG)
#*
#***************************************
GET_LOG_PROPOSAL_Q_UNI_VAR <- function(mcmc_samples, epidemic_data, n_samples, 
                               dof = 3, prob = 0.95, FLAG_PRIORS = list(EXP_PRIOR = TRUE)) {    # FLAG_DIM = list(UNI_VAR = FALSE, MULTI_VAR = FALSE)
  
  'Get proposal q '
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  
  #SAMPLING SIZE
  samp_size_proposal = prob*n_samples; samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #*******
  #THETA SAMPLES: PROPOSAL + PRIOR
  mean_mcmc = mean(mcmc_samples); sd_mcmc = sd(mcmc_samples)
  theta_samples_proposal = sd_mcmc*rt(samp_size_proposal, df = dof) + mean_mcmc 
  
  #PRIORS
  if(FLAG_PRIORS$EXP_PRIOR){
    theta_samples_prior = c(rexp(samp_size_prior))
  }
  
  #MIXTURE Q
  theta_samples = c(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dt((theta_samples - mean_mcmc)/sd_mcmc, df = 1, log = TRUE) - log(sd_mcmc) 
  
  #PRIOR
  if(FLAG_PRIORS$EXP_PRIOR){
    log_prior = dexp(theta_samples, log = TRUE) 
  }
  
  #LOG Q
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_prior)) #CALCULATE WITH LOG SUM EXP TRICK ASWELL & SEE IF MATCH
  
  #LOG SUM EXP TRICK TO GET LOG_Q (MATCH) 
  max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_prior)
  log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_prior - max_el))
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps)
}


#*************************************
#1b. MULTI PARAMETER MODEL - PROPOSAL Q
#*************************************
GET_LOG_PROPOSAL_Q_MULTI_DIM <- function(mcmc_samples, epidemic_data,  
                                         n_samples, n_dims = 3, dof = 3, prob = 0.95, 
                                         FLAG_DIMS = list(dims_two = FALSE, dims_three = TRUE),
                                         FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                             SSIB = FALSE, SSIC = FALSE),
                                         FLAG_PRIORS = list(EXP_PRIOR = TRUE)) {               
  
  #PARAMETERS REQUIRED
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR
  means = colMeans(mcmc_samples)
  #theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) + means 
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal) 
  
  if(FLAGS_MODELS$SSEB){
    theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = 3) 
  } else if (FLAGS_MODELS$SSNB){
    theta_samples_prior = matrix(c(r(samp_size_prior), r(samp_size_prior)))
  }

  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  print(paste0('mean_mcmc = ', means))
  #print(paste0('mean proposal = ', colMeans(theta_samples_proposal)))
  #print(paste0('mean prior = ', colMeans(theta_samples_prior)))
  
  #DEFENSE MIXTURE
  #log_proposal = dmvt(theta_samples - means, sigma = cov(mcmc_samples), df = dof, log = TRUE)
  log_proposal = dmvt(theta_samples - matrix( rep(means, each = n_samples), ncol = n_dims),
                      sigma = cov(mcmc_samples), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  = 2, f(x,y) = 0.2) for exmaples
  log_prior = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) + dexp((theta_samples[,3] - 1), log = TRUE)
  
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_prior)) #1 x n_samples
  
  # print(paste0('mean_mcmc = ', means))
  # print(paste0('mean log_proposal = ', mean(log_proposal)))
  # print(paste0('mean prior = ', mean(log_prior, na.rm = TRUE)))
  # print(paste0('1 log_q = ', log_q))
  # 
  # max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_prior)
  # log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_prior - max_el))
  # log_q_s = LOG_SUM_EXP(log_q2) #LOG SUM EXP OF TWO COMPONENTS - See if the same
  # print(paste0('2 log_q_s = ', log_q_s))
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps)  
}