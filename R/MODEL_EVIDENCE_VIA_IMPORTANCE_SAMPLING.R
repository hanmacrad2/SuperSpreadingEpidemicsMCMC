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
#' @export
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
#' @export
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
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_prior = log_prior)
  
  return(imp_samp_comps)
}


#*************************************
#1b. MULTI PARAMETER MODEL - PROPOSAL Q
#*************************************
#' @export
GET_LOG_PROPOSAL_Q_MULTI_DIM <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                                         n_samples, dof = 3, prob = 0.95, 
                                         priors_ssnb = list(pk_ga_shape = 0.001, pk_ga_rte = 0.001, pr0_unif = c(1.0,4),
                                                       p_prob_unif = c(0,1)),
                                         FLAG_PRIORS = list(SSEB_EXP_PRIOR = TRUE)) {               
  
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
  
  neg_count <- apply(theta_samples_proposal, 2, function(x) sum(x < 0))
  print(paste0('neg_count, ', neg_count))
  
  #PRIORS
  #print(paste0('FLAGS_MODELS$SSNB', FLAGS_MODELS))
  if(FLAGS_MODELS$SSEB & FLAG_PRIORS$SSEB_EXP_PRIOR){
    n_dim = 3 #PUT IN DICTIONARY
    theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = n_dim) 
    
  } else if (FLAGS_MODELS$SSNB){
    n_dim = 2
    theta_samples_prior = matrix(c(rgamma(samp_size_prior, shape = priors_ssnb$pk_ga_shape, rate = priors_ssnb$pk_ga_rte),
                                   runif(samp_size_prior,  min = priors_ssnb$pr0_unif[1], max = priors_ssnb$pr0_unif[2])), ncol = n_dim)
  }
  
  #PROPOSAL SAMPLES
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dmvt(theta_samples - matrix(rep(means, each = n_samples), ncol = n_dim),
                      sigma = cov(mcmc_samples), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #PRIORS
  if(FLAGS_MODELS$SSEB  & FLAG_PRIORS$SSEB_EXP_PRIOR){
    log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE) 
  } else if (FLAGS_MODELS$SSNB){
    log_priors = dgamma(theta_samples[,1], shape = priors_ssnb$pk_ga_shape, rate = priors_ssnb$pk_ga_rte, log = TRUE) +
      dunif(theta_samples[, 2], min = priors_ssnb$pr0_unif[1], max = priors_ssnb$pr0_unif[2], log = TRUE)
  }
 
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_priors)) #1 x n_samples
  
  #max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_priors)
  #log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_priors - max_el))
  #log_q_s = LOG_SUM_EXP(log_q2) #LOG SUM EXP OF TWO COMPONENTS - See if the same
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_priors = log_priors)
  
  return(imp_samp_comps)  
}

#*********************************************************
#*
#* 2. GET P_HATS ESTIMATE OF MODEL EVIDENCE (LOG)
#*
#************************************************************
#' @export
GET_LOG_P_HAT_BASELINE <- function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PROPOSAL, PRIORS
  imp_samp_comps = GET_LOG_PROPOSAL_Q_UNI_VAR(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples 
  log_q = imp_samp_comps$log_q; log_priors = imp_samp_comps$log_priors
  
  #PRIORS 
  log_priors = dexp(theta_samples, log = TRUE) 
  
  #LOG SUM EXP (LOOP)
  vector_log_sum_exp = rep(NA, n_samples)
  
  for(i in 1:n_samples){         
    
    loglike = LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])
    vector_log_sum_exp[i] = loglike + log_priors[i] - log_q[i]
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
  
  return(log_p_hat)
}

#******************************************************************
#* 2b. Estimate of model evidence (P_hat) for SSEB model 
#******************************************************************
#' @export
GET_LOG_P_HAT <- function(mcmc_samples, epidemic_data,
                          FLAGS_MODELS, n_samples = 200,
                          priors_ssnb = list(pk_ga_shape = 0.001, pk_ga_rte = 0.001, pr0_unif = c(1.0,4),
                                             p_prob_unif = c(0,1)),
                          FLAG_PRIORS = list(SSEB_EXP_PRIOR = TRUE)) {    
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(0, n_samples)
  lambda_vec = get_lambda(epidemic_data); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q_MULTI_DIM(mcmc_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_priors = imp_samp_comps$log_priors
  
  #SSEB MODEL 
  if (FLAGS_MODELS$SSEB & FLAG_PRIORS$SSEB_EXP_PRIOR) {
    
    #GET ESTIMATE
    for (i in 1:n_samples) {
      count_samps = 0
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3])
      
      vector_estimate_terms = loglike + log_priors[i] - log_q[i]
      
      if(!is.na(vector_estimate_terms) ||  !is.infinite(vector_estimate_terms)){
        
        vector_estimate_terms[i] =  vector_estimate_terms
        #count_samps = count_samps + 1
        #print(paste0('count_samps: ', count_samps))
        
      }
    }
    #print(paste0('total count_samps: ', count_samps))
    
    #SSNB MODEL
  } else if (FLAGS_MODELS$SSNB) {
    
    for (i in 1:n_samples) {
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_SSNB(epidemic_data, lambda_vec, theta_samples[i,]) 
      
      vector_estimate_terms[i] = loglike + log_priors[i] - log_q[i]
    }
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  log_p_hat2 = -log(n_samples) + log(sum(exp(vector_estimate_terms)))
  print(paste0('log_p_hat2 = ', log_p_hat2))
  
  return(log_p_hat)
}


#******************************************************************************
#*
# 3. LOAD MCMC + GET P_HAT ESTIMATES MODEL EVIDENCE ESTIMATES
#*
#******************************************************************************
#' @export
LOAD_MCMC_GET_P_HAT <- function(epidemic_data, OUTER_FOLDER, run = run, n_repeats = n_repeats,  BURN_IN = FALSE, burn_in_pc = 0.2,
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                    SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Load mcmc samples 2. Get estimate'
  
  #Parameters
  estimates_vec = rep(NA, n_repeats) 
  
  if (FLAGS_MODELS$BASE){
    
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      
      #READ SAMPLES
      print(paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      print(mcmc_output$r0_vec)
      #print(paste0(CURRENT_FOLDER, 'mcmc_base_', i))
      #mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_base_', i))
     
      
      #LOG LIKE (BURN-IN):
      if(BURN_IN){
        mcmc_output$r0_vec = mcmc_output$r0_vec[(burn_in_pc*length(mcmc_output$r0_vec)):length(mcmc_output$r0_vec)]
        print(length(mcmc_output$r0_vec))
      }
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      log_phat = GET_LOG_P_HAT_BASELINE(mcmc_output$r0_vec, epidemic_data) 
      estimates_vec[i] = log_phat
      print(estimates_vec)
    }
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/phat_ests_base_', run, '.rds'))
    
  } else if(FLAGS_MODELS$SSEB){
    
    model_type = 'sseb'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      #mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i))
      mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data, FLAGS_MODELS = FLAGS_MODELS)
      
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/phat_ests_', model_type, '_', run, '.rds'))
    
  } else if (FLAGS_MODELS$SSNB){
    
    model_type = 'ssnb'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', i ,'.rds'))
      
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_samples =  mcmc_output$ssnb_params_matrix 
      print(mcmc_samples[1])
      #LOG LIKE (BURN-IN):
      if(BURN_IN){
        mcmc_length = dim(mcmc_output$ssnb_params_matrix)[1]
        mcmc_samples = mcmc_output$ssnb_params_matrix[(burn_in_pc*mcmc_length + 1):mcmc_length,]
        print(dim(mcmc_samples))
      }
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data, FLAGS_MODELS)
      estimates_vec[i] = phat_estimate                        
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, '/phat_ests_', model_type, '_', run, '.rds' ))
    
  }
  
  return(estimates_vec) 
}

#***************************** 
# LOAD MCMC VERSION 2
#' @export
LOAD_MCMC_GET_P_HAT_II <- function(matrix_data, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                    SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Load mcmc samples 2. Get estimate'
  
  #Parameters
  estimates_vec = rep(NA, n_repeats) 
  
  if (FLAGS_MODELS$BASE){
    
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/')
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      
      #READ SAMPLES
      epidemic_data = matrix_data[i, ]
      print(epidemic_data)
    
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      print(mean(mcmc_output$r0_vec))
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      log_phat = GET_LOG_P_HAT_BASELINE(mcmc_output$r0_vec, epidemic_data) 
      estimates_vec[i] = log_phat
      print(estimates_vec)
    }
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, 'phat_ests_base_', run, '.rds'))
    
  } else if(FLAGS_MODELS$SSEB){
    
    model_type = 'sseb'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/')
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      
      #MCMC
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      if( sum(is.na(mcmc_output$log_like_vec)) == 0){
        
        print(paste0('working i = ', i))
        mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)
        #DATA
        epidemic_data = matrix_data[i, ]
        
        #GET PHAT ESTIMATE OF MODEL EVIDENCE
        phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data, FLAGS_MODELS = FLAGS_MODELS)
        
        estimates_vec[i] = phat_estimate
        print(estimates_vec)
      }

      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, 'phat_ests_', model_type, '_', run, '.rds'))
    
  } else if (FLAGS_MODELS$SSNB){
    
    model_type = 'ssnb'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type),  '/')
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      epidemic_data = matrix_data[i, ]
      mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
      mcmc_samples =  mcmc_output$ssnb_params_matrix 
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data, FLAGS_MODELS)
      estimates_vec[i] = phat_estimate                        
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_FOLDER, 'phat_ests_', model_type, '_', run, '.rds' ))
    
  }
  
  return(estimates_vec) 
}
