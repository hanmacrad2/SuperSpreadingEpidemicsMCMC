#********************************************************************
#MODEL EVIDENCE SSI MODEL (DATA AUG)

#************************
# 1. MODEL EVIDENCE (P HAT)
#************************
GET_LOG_MODEL_EVIDENCE_SSI <- function(mcmc_output, EPI_DATA, 
                                       FLAGS_MODELS = list(BASE = FALSE, SSE = FALSE,
                                                           SSI = TRUE, SSEB = FALSE, SSIB = FALSE), n_samples = 10000){
  
  'Estimate of model evidence for SSI model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(NA, n_samples)
  
  #vector_estimate_terms = c()
  lambda_vec = get_lambda(EPI_DATA); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q_SSI(mcmc_output, EPI_DATA, FLAGS_MODELS, n_samples) 
  
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  infectivity_vec = GET_INFECT_GAMMA_CURVE(EPI_DATA) 
  num_etas = length(EPI_DATA)-1 
  count_ok = 0
  
  for (i in 1:n_samples) {
    
    if (log_prior_density[i] > -Inf) {
      
      loglike = LOG_LIKE_SSI(EPI_DATA, infectivity_vec, theta_samples[i, 1:2],
                             theta_samples[i, 3:dim(theta_samples)[2]], FLAG_MCMC = FALSE) 
    } else {
      loglike = -Inf
    }
    
    vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
    
  }
  
  #ESTIMATE OVER SUM
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  return(log_p_hat)
  
  
}

#********************************************************************
#*
# GET PROPOSAL DENSITY SSI MODEL
#*
#********************************************************************
GET_LOG_PROPOSAL_Q_SSI <- function(mcmc_output, EPI_DATA, FLAGS_MODELS,
                                    n_samples, dof = 3, prob = 0.5) { #prob = 0.95 0.9999
  
  #PARAMETERS REQUIRED 
  lambda_vec = get_lambda(EPI_DATA)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; 
  samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #MCMC SAMPLES
  mcmc_samples = cbind(mcmc_output$ssi_params_matrix, mcmc_output$eta_matrix) #eta_nonzero_cols)  
  n_dim = dim(mcmc_samples)[2] 
  means = colMeans(mcmc_samples[,])

  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal)
  
  theta_samples_prior = GET_PRIOR_IMPORTANCE_SAMPLES(EPI_DATA, samp_size_prior, FLAGS_MODELS)
  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  #DENSITY OF PROPOSAL
  wh_non_zero = which(EPI_DATA[1:(length(EPI_DATA)-1)]!= 0)
  
  log_proposal_density = dmvt(theta_samples[,2+wh_non_zero,drop=FALSE] - matrix_means[,2+wh_non_zero,drop=FALSE], #wh_non_zero: Include wh here to only include non zero eta columns 
                              sigma = cov(mcmc_samples[,2+wh_non_zero,drop=FALSE]), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples

  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples, EPI_DATA, FLAGS_MODELS) #n_dim
  
  #PROPOSAL 
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_prior_density)) #1 x n_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior_density = log_prior_density)
  
  return(imp_samp_comps)  
}


#************************
#1. LOAD MCMC FOR MODEL EVIDENCE 
#************************
LOAD_MCMC_GET_SSI_MODEL_EV <- function(EPI_DATA, OUTER_FOLDER, 
                                        run = run, n_repeats = n_repeats, 
                                        start = 1){
  
  #Parameters
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  estimates_vec = rep(NA, n_repeats) 
  
  model_type = 'SSI'; #print(model_type)
  CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  print(CURRENT_FOLDER)
  
  for (i in start:n_repeats){ 
    mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
    
    #GET PHAT ESTIMATE OF MODEL EVIDENCE
    phat_estimate = GET_LOG_MODEL_EVIDENCE_SSI(mcmc_output, EPI_DATA) 
    estimates_vec[i] = phat_estimate                        
    print(estimates_vec)
  }
  
  #SAVE ESTIMATES
  saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
  
  return(estimates_vec) 
}
