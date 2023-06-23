#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)
#library(mvtnorm)

#*************************************
#* FUNCTIONS 
#* ***********************************
LOG_SUM_EXP <- function(vectorX){
  
  #REMOVE NA VALUES
  vectorX = na.omit(vectorX)
  
  max_val = max(vectorX)
  
  out = max_val + log(sum(exp(vectorX - max_val)))
  
  return(out)
}

#********************************************************************
#*
#1. GET PROPOSAL DENSITY 
#*
#********************************************************************
GET_LOG_PROPOSAL_Q <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                                         n_samples, dof = 3, prob = 0.95) {
  
  #PARAMETERS REQUIRED 
  n_dim = dim(mcmc_samples)[2] 
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
  theta_samples_prior = GET_PRIOR_THETA_SAMPLES(epidemic_data, samp_size_prior, n_dim, FLAGS_MODELS)
  
  print(paste0('theta samps proposal dim: ', dim(theta_samples_proposal)))
  print(paste0('theta samps prior dim: ', dim(theta_samples_prior)))
  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  #DENSITY OF PROPOSAL
  log_proposal_density = dmvt(theta_samples - matrix_means,
                      sigma = cov(mcmc_samples), df = dof) #log = TRUE log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples, epidemic_data,
                                            samp_size_prior, n_dim, FLAGS_MODELS)
  
  #PROPOSAL 
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_prior_density)) #1 x n_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior_density = log_prior_density)
  
  return(imp_samp_comps)  
}


#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE (LOG) (P_HAT)
#*
#*****************************************************************************************
GET_LOG_MODEL_EVIDENCE <- function(mcmc_samples, epidemic_data,
                          FLAGS_MODELS, n_samples = 1000) {   
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(NA, n_samples)
  lambda_vec = get_lambda(epidemic_data); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #SSEB MODEL 
  if (FLAGS_MODELS$SSEB){
    
    #GET ESTIMATE
    for (i in 1:n_samples) {
      
      loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3])
      
      vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
    }
    
    #SSNB MODEL
  } else if (FLAGS_MODELS$SSNB) {
    
    for (i in 1:n_samples) {
      
      loglike = LOG_LIKE_SSNB(epidemic_data, lambda_vec, theta_samples[i,]) 
      
      vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
    }
  } else if (FLAGS_MODELS$SSIR) {
    
    infectivity_vec = get_infectious_curve(epidemic_data)
    num_etas = length(epidemic_data)-1
    
    for (i in 1:n_samples) {
      loglike = LOG_LIKE_SSIR(epidemic_data, infectivity_vec, theta_samples[i, 1:2],  theta_samples[i, 3:(2+num_etas)]) 
      
      vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
    }
    
  } else if (FLAGS_MODELS$SSIB){
    
    #GET ESTIMATE
    for (i in 1:n_samples) {
      
      if (log_prior_density[i] > -Inf){

        loglike = LOG_LIKE_SSIB(epidemic_data, theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3])
      
        vector_estimate_terms[i] = loglike + log_prior_density[i] - log_q[i]
      } else {
        vector_estimate_terms[i] = -Inf
      }
    }
    
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  return(log_p_hat)
}

#******************************************************************************
#*
# 3. LOAD MCMC + GET P_HAT ESTIMATES MODEL EVIDENCE ESTIMATES
#*
#******************************************************************************
LOAD_MCMC_GET_MODEL_EVIDENCE <- function(epidemic_data, OUTER_FOLDER, 
                                         run = run, n_repeats = n_repeats, start = 1, 
                                         beta_ssib = 1,
                                         PRIORS_LIST = list(),
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                    SSIB = FALSE, SSIR = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Load mcmc samples 2. Get estimate'
  
  #Parameters
  estimates_vec = rep(NA, n_repeats) 
  #estimates_vec = rep(NA, n_repeats - start) 
  
  if (FLAGS_MODELS$BASE){
    
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in start:n_repeats){
      
      print(paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))

      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      log_phat = GET_LOG_MODEL_EVIDENCE_BASELINE(mcmc_output$r0_vec, epidemic_data, PRIORS_LIST) 
      estimates_vec[i] = log_phat
      print(estimates_vec)
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
    
  } else if(FLAGS_MODELS$SSEB){
    
    model_type = 'sseb'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in start:n_repeats){
      
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data, FLAGS_MODELS = FLAGS_MODELS)
      
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
    
  } else if (FLAGS_MODELS$SSNB){
    
    model_type = 'ssnb'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in start:n_repeats){
      
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_samples =  mcmc_output$ssnb_params_matrix 
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data, FLAGS_MODELS)
      estimates_vec[i] = phat_estimate                        
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
    
  } else if (FLAGS_MODELS$SSIR) {
    
    model_type = 'ssir'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    print(CURRENT_FOLDER)
    
    for (i in start:n_repeats){ #start:n_repeats
      
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_samples = cbind(mcmc_output$ssir_params_matrix, mcmc_output$eta_matrix)
      print(paste0('dim of mcmc', dim(mcmc_samples)))

      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data, FLAGS_MODELS)
      estimates_vec[i] = phat_estimate                        
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
    
  } else if (FLAGS_MODELS$SSIB){
    
    model_type = 'ssib'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    print(CURRENT_FOLDER)
    
    for (i in start:n_repeats){ #start:n_repeats

      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      #mcmc_samples =  matrix(c(mcmc_output[["a_vec"]], mcmc_output[["b_vec"]], mcmc_output[["c_vec"]]), ncol = 3)
      #mcmc_samples =  matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)

      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_output, epidemic_data, beta = beta_ssib, FLAGS_MODELS)
      estimates_vec[i] = phat_estimate                        
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
    
  }
  
  return(estimates_vec) 
}
