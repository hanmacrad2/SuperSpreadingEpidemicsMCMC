#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)

#********************************************************************
#*
#1. GET IMPORTANCE SAMPLING PROPOSAL (LOG)
#*
#********************************************************************

GET_LOG_PROPOSAL_Q_UNI_VAR <- function(mcmc_samples, epidemic_data, n_samples, 
                                       PRIORS_LIST = list(EXP_PRIOR = FALSE, GAMMA_PRIOR = TRUE), 
                                       dof = 3, prob = 0.95, r0_sim = 1.6) {    # FLAG_DIM = list(UNI_VAR = FALSE, MULTI_VAR = FALSE)
  
  'Get proposal q '
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  
  #SAMPLING SIZE
  samp_size_proposal = prob*n_samples; samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  print(paste0('samp_size_proposal = ', samp_size_proposal))
  
  #*******
  #THETA SAMPLES: PROPOSAL + PRIOR
  mean_mcmc = mean(mcmc_samples); sd_mcmc = sd(mcmc_samples)
  theta_samples_proposal = sd_mcmc*rt(samp_size_proposal, df = dof) + mean_mcmc 
  
  #PRIORS
  if(PRIORS_LIST$EXP_PRIOR){
    theta_samples_prior = c(rexp(samp_size_prior))
  } else if (PRIORS_LIST$GAMMA_PRIOR) {
    gamma_shape = 2 
    theta_samples_prior = c(rgamma(samp_size_prior, shape = 2, scale = gamma_shape*r0_sim))
  }
  
  #MIXTURE Q
  theta_samples = c(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dt((theta_samples - mean_mcmc)/sd_mcmc, df = 1, log = TRUE) - log(sd_mcmc) 
  
  #PRIOR
  if(PRIORS_LIST$EXP_PRIOR){
    log_prior = dexp(theta_samples, log = TRUE) 
  } else if (PRIORS_LIST$GAMMA_PRIOR){
    gamma_shape = 2
    log_prior = dgamma(theta_samples, shape = gamma_shape, scale = gamma_shape*r0_sim, log = TRUE)
  }
  
  #LOG Q
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_prior)) #CALCULATE WITH LOG SUM EXP TRICK ASWELL & SEE IF MATCH
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior = log_prior)
  
  return(imp_samp_comps)
}

#**************************************************************************************
#*
#* 2. GET P_HATS ESTIMATE OF MODEL EVIDENCE (LOG)
#*
#*****************************************************************************************
GET_LOG_MODEL_EVIDENCE_BASELINE <- function(mcmc_samples, epidemic_data, PRIORS_LIST, 
                                            n_samples = 1000, r0_sim = 1.6) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PROPOSAL, PRIORS
  imp_samp_comps = GET_LOG_PROPOSAL_Q_UNI_VAR(mcmc_samples, epidemic_data, n_samples, PRIORS_LIST)
  theta_samples = imp_samp_comps$theta_samples 
  log_q = imp_samp_comps$log_q; log_priors = imp_samp_comps$log_priors
  
  #PRIORS 
  if (PRIORS_LIST$GAMMA_PRIOR){
    gamma_shape = 2
    log_priors = dgamma(theta_samples,  shape = gamma_shape, scale = gamma_shape*r0_sim, log = TRUE)
    
  } else {
    log_priors = dexp(theta_samples, log = TRUE)
  }
  
  #LOG SUM EXP (LOOP)
  vector_log_sum_exp = rep(NA, n_samples)
  
  for(i in 1:n_samples){         
    
    loglike = LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])
    vector_log_sum_exp[i] = loglike + log_priors[i] - log_q[i]
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
  
  return(log_p_hat)
}
