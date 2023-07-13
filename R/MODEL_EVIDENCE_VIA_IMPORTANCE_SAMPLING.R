#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)

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
                                         n_samples, dof = 3, prob = 0.9999) { #prob = 0.95 0.9999
  
  #PARAMETERS REQUIRED 
  n_dim = dim(mcmc_samples)[2] 
  print(paste0('n_dim:', n_dim))
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
  
  #print(paste0('dim theta samps proposal dim: ', dim(theta_samples_proposal)))
  #print(paste0('dim theta samps prior dim: ', dim(theta_samples_prior)))
  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  #print(paste0('dim theta_samples', dim(theta_samples)))
  #print(paste0('dim matrix_means', dim(matrix_means)))
  
  #DENSITY OF PROPOSAL
  log_proposal_density = dmvt(theta_samples - matrix_means,
                      sigma = cov(mcmc_samples), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
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


#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE (LOG) (P_HAT)
#
#*****************************************************************************************


#******************************************************************************
#*
# 3. LOAD MCMC + GET P_HAT ESTIMATES MODEL EVIDENCE ESTIMATES
#*
#******************************************************************************
LOAD_MCMC_GET_MODEL_EVIDENCE <- function(epidemic_data, OUTER_FOLDER, 
                                         run = run, n_repeats = n_repeats, 
                                         start = 1, beta_ssib = 1000,
                                         num_is_samps = 10000,
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                    SSIB = FALSE, SSIR = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Load mcmc samples 2. Get estimate'
  
  #Parameters
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  estimates_vec = rep(NA, n_repeats) 
  #estimates_vec = rep(NA, n_repeats - start) 
  
  if (FLAGS_MODELS$BASE){
    
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in start:n_repeats){
      
      print(paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))

      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      log_phat = GET_LOG_MODEL_EVIDENCE_BASELINE(mcmc_output$r0_vec, epidemic_data) 
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
      
      #ETA ONLY IF > 0
      eta_nonzero <- apply(mcmc_output$eta_matrix != 0, 2, any)
      eta_nonzero_cols = mcmc_output$eta_matrix[, eta_nonzero]
      
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
      phat_estimate = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_output, epidemic_data, beta = beta_ssib,
                                                  num_is_samps = num_is_samps, FLAGS_MODELS)
      estimates_vec[i] = phat_estimate                        
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
    
  }
  
  return(estimates_vec) 
}
