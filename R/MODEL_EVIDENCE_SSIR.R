#********************************************************************
#MODEL EVIDENCE SSIR MODEL (DATA AUG)
#********************************************************************


#********************************************************************
#*
#1. GET PROPOSAL DENSITY SSIR MODEL
#*
#********************************************************************
GET_LOG_PROPOSAL_Q_SSIR <- function(mcmc_output, EPI_DATA, FLAGS_MODELS,
                                    n_samples, dof = 3, prob = 0.5) { #prob = 0.95 0.9999
  
  #browser()
  #PARAMETERS REQUIRED 
  lambda_vec = get_lambda(EPI_DATA)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; 
  samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #ETA REMOVE 0
  #wh_nonzero <- apply(mcmc_output$eta_matrix != 0, 2, any)
  #eta_nonzero_cols = mcmc_output$eta_matrix[, wh_nonzero]
  mcmc_samples = cbind(mcmc_output$ssir_params_matrix, mcmc_output$eta_matrix) #eta_nonzero_cols)  
  n_dim = dim(mcmc_samples)[2] 
  
  #THETA SAMPLES: PROPOSAL + PRIOR (FROM PARAMETRIC APPROXIMATION)
  means = colMeans(mcmc_samples[,])

  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal)
  
  theta_samples_prior = GET_PRIOR_THETA_SAMPLES(EPI_DATA, samp_size_prior, n_dim, FLAGS_MODELS)
  
  #which_small = apply(theta_samples_prior < 0.005, 2, any)
  #print(theta_samples_prior[,which_small])
  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  #DENSITY OF PROPOSAL ** ADDED!! 14/07/23
  wh_non_zero = which(EPI_DATA[1:(length(EPI_DATA)-1)]!= 0)
  
  log_proposal_density = dmvt(theta_samples[,2+wh_non_zero,drop=FALSE] - matrix_means[,2+wh_non_zero,drop=FALSE], #wh_non_zero: Include wh here to only include non zero eta columns 
                              sigma = cov(mcmc_samples[,2+wh_non_zero,drop=FALSE]), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples

  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples, EPI_DATA,
                                            samp_size_prior, n_dim, FLAGS_MODELS)
  
  #PROPOSAL 
  
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_prior_density)) #1 x n_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior_density = log_prior_density)
  
  return(imp_samp_comps)  
}

#************************
# 2. MODEL EVIDENCE (P HAT)
#************************
GET_LOG_MODEL_EVIDENCE_SSIR <- function(mcmc_output, EPI_DATA, 
                                        FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                            SSIB = FALSE, SSIR = TRUE), n_samples = 10000){
  
  #browser()
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  #eta_nonzero <- apply(mcmc_output$eta_matrix != 0, 2, any)
  #eta_nonzero_cols = mcmc_output$eta_matrix[, eta_nonzero]
  #mcmc_samples = cbind(mcmc_output$ssir_params_matrix, mcmc_output$eta_matrix) #eta_nonzero_cols)
  
  vector_estimate_terms = rep(NA, n_samples)
  
  #vector_estimate_terms = c()
  lambda_vec = get_lambda(EPI_DATA); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q_SSIR(mcmc_output, EPI_DATA, FLAGS_MODELS, n_samples) 
  
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  #print(paste0('dim theta samples', dim(theta_samples)))
  
  infectivity_vec = GET_INFECT_GAMMA_CURVE(EPI_DATA) 
  num_etas = length(EPI_DATA)-1 # dim(theta_samples)[2] - 2
  count_ok = 0
  
  for (i in 1:n_samples) {
    
    if (log_prior_density[i] > -Inf) {
    
      loglike = LOG_LIKE_SSIR(EPI_DATA, infectivity_vec, theta_samples[i, 1:2],
                              theta_samples[i, 3:dim(theta_samples)[2]], FLAG_MCMC = FALSE) 
      count_ok = count_ok + 1
      #browser()
    } else {
      loglike = -Inf
    }
    
    #browser()
    vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
    
   # if (loglike > 0){
      #print(paste0('vector_estimate_terms[i]', vector_estimate_terms[i]))
      #browser()
  #  }
    #print(paste0('vector_estimate_terms[i]', vector_estimate_terms[i]))
    
    
  }
  
  
  #ESTIMATE OVER SUM
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  # if(log_p_hat > 1){
  #   print(paste0('POSTIVE log_p_hat = ', log_p_hat))
  #   print(paste0('loglike = ', loglike))
  #   print(paste0('vec terms', vector_estimate_terms))
  # } 
  
  return(log_p_hat)
  
  
}


#************************
#1. LOAD MCMC FOR MODEL EVIDENCE 
#************************
LOAD_MCMC_GET_SSIR_MODEL_EV <- function(EPI_DATA, OUTER_FOLDER, 
                                        run = run, n_repeats = n_repeats, 
                                        start = 1){
  
  #Parameters
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  estimates_vec = rep(NA, n_repeats) 
  
  model_type = 'ssir'; print(model_type)
  CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  print(CURRENT_FOLDER)
  
  for (i in start:n_repeats){ 
    #print(paste0('i = ', i))
    mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
    
    #GET PHAT ESTIMATE OF MODEL EVIDENCE
    phat_estimate = GET_LOG_MODEL_EVIDENCE_SSIR(mcmc_output, EPI_DATA) 
    estimates_vec[i] = phat_estimate                        
    print(estimates_vec)
  }
  
  #SAVE ESTIMATES
  saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
  
  return(estimates_vec) 
}

#MULITPLE RESULTS
GET_MAT_RESULTS_SSIR <- function(mat_ssir, n_reps = 10){
  
  for (i in 3:n_reps) {
    mat_ssir[,i] = LOAD_MCMC_GET_SSIR_MODEL_EV(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats)
  }
  
  return(mat_ssir)
}


#RUN
#PARAMS
#run = 1; n_repeats = 5
# model_ev_ssir6 = LOAD_MCMC_GET_SSIR_MODEL_EV(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats)
# mean(model_ev_ssir6)
# sd(model_ev_ssir6)
