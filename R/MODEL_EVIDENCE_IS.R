#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#*************************************
#* FUNCTIONS 
#* ***********************************
LOG_SUM_EXP <- function(vectorX){
  
  #REMOVE NA VALUES
  vectorX = na.omit(vectorX)
  vectorX = vectorX[!is.nan(vectorX)]
  vectorX = vectorX[!is.infinite(vectorX)]
  
  max_val = max(vectorX)
  
  out = max_val + log(sum(exp(vectorX - max_val)))
  
  # if (is.nan(out)){
  #   print(paste0('max val', max_val))
  # }
  return(out)
}

#********************************************************************
#*
#1. GET PROPOSAL DENSITY 
#*
#********************************************************************
GET_LOG_PROPOSAL_Q <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                               n_samples, dof = 3, prob = 0.95) { #0.5  #prob = 0.95 0.9999
  
  #PARAMETERS REQUIRED 
  n_dim = dim(mcmc_samples)[2] 
  #print(paste0('n_dim:', n_dim))
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; 
  samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR (FROM PARAMETRIC APPROXIMATION)
  means = colMeans(mcmc_samples)
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal) 
  theta_samples_prior = GET_PRIOR_IMPORTANCE_SAMPLES(epidemic_data, samp_size_prior, FLAGS_MODELS)
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  #DENSITY OF PROPOSAL
  log_proposal_density = dmvt(theta_samples - matrix_means,
                              sigma = cov(mcmc_samples), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples, epidemic_data, FLAGS_MODELS)
  
  #PROPOSAL 
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_prior_density)) #1 x n_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior_density = log_prior_density)
  
  return(imp_samp_comps)  
} #0.001

#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE (LOG) (P_HAT)
#
#*****************************************************************************************
GET_LOG_MODEL_EVIDENCE <- function(mcmc_samples, epidemic_data,
                                   FLAGS_MODELS = list(BASELINE = FALSE, SSE = FALSE,
                                                       SSI = FALSE, SSEB = FALSE, SSIB = FALSE), n_samples = 10000) {   
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(NA, n_samples)
  lambda_vec = get_lambda(epidemic_data);
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #GET ESTIMATE
  for (i in 1:n_samples) {
    
    if (log_prior_density[i] > -Inf ) {
      
      #SSEB MODEL 
      if (FLAGS_MODELS$SSEB){
        loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],
                                theta_samples[i, 2], theta_samples[i, 3])
        
        #SSE MODEL 
      } else if (FLAGS_MODELS$SSE) {
        
        loglike = LOG_LIKE_SSE(epidemic_data, lambda_vec, theta_samples[i,]) 
        
      } 
      
    } else {
      
      loglike = 0
    }
    
    vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
    
  }
  
  #ESTIMATE OVER SUM
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  return(log_p_hat)
  
} 