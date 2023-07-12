#MODEL EVIDENCE SSIR MODEL (DATA AUG)

#************************
#1. LOAD MCMC FOR MODEL EVIDENCE 
#************************
LOAD_MCMC_GET_SSIR_MODEL_EV <- function(EPI_DATA, OUTER_FOLDER, 
                                        run = run, n_repeats = n_repeats, start = 1, num_is_samps = 10000,
                                        FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                            SSIB = FALSE, SSIR = truehist())){
  
  
  model_type = 'ssir'; print(model_type)
  CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  print(CURRENT_FOLDER)
  
  for (i in start:n_repeats){ #start:n_repeats
    
    mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
    
    # #ETA ONLY IF > 0
    # eta_nonzero <- apply(mcmc_output$eta_matrix != 0, 2, any)
    # eta_nonzero_cols = mcmc_output$eta_matrix[, eta_nonzero]
    
    mcmc_samples = cbind(mcmc_output$ssir_params_matrix, mcmc_output$eta_matrix)
    print(paste0('dim of mcmc', dim(mcmc_samples)))
    
    #GET PHAT ESTIMATE OF MODEL EVIDENCE
    phat_estimate = GET_LOG_MODEL_EVIDENCE_SSIR(mcmc_output, EPI_DATA, FLAGS_MODELS)
    estimates_vec[i] = phat_estimate                        
    print(estimates_vec)
    
  }
  
  #SAVE ESTIMATES
  saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
  
}

#************************
# 2. MODEL EVIDENCE (P HAT)
#************************
GET_LOG_MODEL_EVIDENCE_SSIR <- function(mcmc_output, EPI_DATA, n_samples = 10000){
  
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  eta_nonzero <- apply(mcmc_output$eta_matrix != 0, 2, any)
  eta_nonzero_cols = mcmc_output$eta_matrix[, eta_nonzero]
  
  mcmc_samples = cbind(mcmc_output$ssir_params_matrix, eta_nonzero_cols)
  
  #vector_estimate_terms = rep(NA, n_samples)
  vector_estimate_terms = c()
  lambda_vec = get_lambda(epidemic_data); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  infectivity_vec = GET_INFECT_GAMMA_CURVE(EPI_DATA) 
  num_etas = dim(theta_samples)[2] - 2 #length(EPI_DATA)-1
  count_inf = 0
  count_nan = 0
  
  for (i in 1:n_samples) {
    
    if (log_prior_density[i] > -Inf && !is.nan(log_prior_density[i])) {
      
      loglike = LOG_LIKE_SSIR(EPI_DATA, infectivity_vec, theta_samples[i, 1:2],
                              theta_samples[i, 3:dim(theta_samples)[2]]) 
      
    } else {
      loglike = 0
    }
    
    vector_estimate_terms_i = loglike + log_prior_density[i] - log_q[i]
    
    vector_estimate_terms = c(vector_estimate_terms, vector_estimate_terms_i)
    
    if (!is.nan(vector_estimate_terms_i) && !is.infinite(vector_estimate_terms_i)) {
      vector_estimate_terms = c(vector_estimate_terms, vector_estimate_terms_i)
    }
  }
}

#********************************************************************
#*
#1. GET PROPOSAL DENSITY SSIR MODEL
#*
#********************************************************************

GET_LOG_PROPOSAL_Q_SSIR <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                               n_samples, dof = 3, prob = 0.95) { #prob = 0.95 0.9999
  
  #PARAMETERS REQUIRED 
  n_dim = dim(mcmc_samples)[2] 
  print(paste0('n_dim:', n_dim))
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; 
  samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #ETA REMOVE 0
  
  #THETA SAMPLES: PROPOSAL + PRIOR (FROM PARAMETRIC APPROXIMATION)
  means = colMeans(mcmc_samples)
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal) 
  theta_samples_prior = GET_PRIOR_THETA_SAMPLES(epidemic_data, samp_size_prior, n_dim, FLAGS_MODELS)
  
  #print(paste0('dim theta samps proposal dim: ', dim(theta_samples_proposal)))
  #print(paste0('dim theta samps prior dim: ', dim(theta_samples_prior)))
  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  #print(paste0('dim theta_samples', dim(theta_samples)))
  #print(paste0('dim matrix_means', dim(matrix_means)))
  
  #DENSITY OF PROPOSAL
  log_proposal_density = dmvt(theta_samples - matrix_means,
                              sigma = cov(mcmc_samples), df = dof) #, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #a = theta_samples[,1:4] - matrix_means[,1:4]
  #b = dmvt(a, sigma = cov(mcmc_samples[,1:4]), df = dof)
  
  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples, epidemic_data,
                                            samp_size_prior, n_dim, FLAGS_MODELS)
  
  #PROPOSAL 
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_prior_density)) #1 x n_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior_density = log_prior_density)
  
  return(imp_samp_comps)  
} #0.001