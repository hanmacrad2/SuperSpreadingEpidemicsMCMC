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
                                       dof = 3, prob = 0.01, r0_sim = 2.0) {    # FLAG_DIM = list(UNI_VAR = FALSE, MULTI_VAR = FALSE)
  
  list_priors = SET_PRIORS()$list_priors
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  
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
  if(PRIORS_USED$BASELINE$r0$EXP){
    theta_samples_prior = c(rexp(samp_size_prior))
  } else if (PRIORS_USED$BASELINE$r0$GAMMA) {
    gamma_shape = 2 
    theta_samples_prior = c(rgamma(samp_size_prior, shape = 2, scale = gamma_shape*r0_sim))
  }
  
  #MIXTURE Q
  theta_samples = c(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dt((theta_samples - mean_mcmc)/sd_mcmc, df = dof, log = TRUE) - log(sd_mcmc) 
  
  #PRIOR
  if(PRIORS_USED$BASELINE$r0$EXP){
    log_prior_density = dexp(theta_samples, log = TRUE) 
  } else if (PRIORS_USED$BASELINE$r0$GAMMA){
    gamma_shape = 2
    log_prior_density = dgamma(theta_samples, shape = gamma_shape, scale = gamma_shape*r0_sim, log = TRUE)
  }
  
  #LOG Q
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_prior_density)) #CALCULATE WITH LOG SUM EXP TRICK ASWELL & SEE IF MATCH
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior_density = log_prior_density)
  
  return(imp_samp_comps)
}

#**************************************************************************************
#*
#* 2. GET P_HATS ESTIMATE OF MODEL EVIDENCE (LOG)
#*
#*****************************************************************************************
GET_LOG_MODEL_EVIDENCE_BASELINE <- function(mcmc_samples, epidemic_data, 
                                            n_samples = 10000, r0_sim = 1.6) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PROPOSAL, PRIORS
  list_priors = SET_PRIORS()$list_priors
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  
  imp_samp_comps = GET_LOG_PROPOSAL_Q_UNI_VAR(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples 
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
 
  #ESS
  sum_weights = 0
  vec_weights_squared = rep(0, times = n_samples)
  
  #PRIORS 
  if (PRIORS_USED$BASELINE$r0$EXP){
    log_prior_density = dexp(theta_samples, log = TRUE)
    
  } else if (PRIORS_USED$BASELINE$r0$GAMMA) {
    gamma_shape = 2
    log_prior_density = dgamma(theta_samples, shape = gamma_shape, scale = gamma_shape*r0_sim, log = TRUE)
    
  }
  
  #LOG SUM EXP (LOOP)
  vector_log_sum_exp = rep(NA, n_samples)
  
  for(i in 1:n_samples){         
    
    #browser()
    if (log_prior_density[i] > -Inf ) {
    loglike = LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])
    
    } else {
      loglike = 0
    }
    vector_log_sum_exp[i] = loglike + log_prior_density[i] - log_q[i] #loglike tiny, log_q a bigger negative and vec_i ends up positive sometines 
    
    sum_weights = sum_weights + exp(vector_log_sum_exp[i])
    vec_weights_squared[i] =  (exp(vector_log_sum_exp[i]))^2
    
    # if ( abs(log_prior_density[i] + log_q[i]) > log(2)){
    #   #print(paste0('vector_estimate_terms[i]', vector_estimate_terms[i]))
    #   browser()
    # }
  }
  # print('CHECK')
  # if (LOG_SUM_EXP(vector_log_sum_exp) > log(n_samples)){
  #   browser() 
  # }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
  
  #ESS = (sum_weights^2)/sum(vec_weights_squared)
  #print(paste0('ESS: ', ESS))
  
  return(log_p_hat)
}
