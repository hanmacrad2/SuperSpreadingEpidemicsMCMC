#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
library(compositions)


#*************************************
#* FUNCTIONS: OTHER 
#* ***********************************
LOG_SUM_EXP <- function(vectorX){
  
  #REMOVE NA VALUES
  vectorX = na.omit(vectorX)
  
  max_val = max(vectorX)
  
  out = max_val + log(sum(exp(vectorX - max_val)))
  
  return(out)
}
#************************
#FUNCTIONS FOR DATA AUGMENTATION MODELS 
#************************
GET_SAMPLES_ETA_PRIORS <- function(param_priors, epidemic_data, samp_size_prior){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  
  time_length =  length(epidemic_data) - 1
  eta_priors_matrix = matrix(nrow = samp_size_prior, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size_prior){
    R0X = param_priors[i, 1]; k = mcmc_samples[i, 2]
    eta_priors_matrix[i, ] = rgamma(time_length, shape = epidemic_data*k, scale = R0X*k)
  }
  
  return(eta_priors_matrix)
}

#************************
#DENSITY
#************************
GET_DENSITY_ETA_PRIORS <- function(theta_samples, epidemic_data){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  
  num_etas =  length(epidemic_data) - 1
  samp_size = dim(theta_samples)[1]
  dim_cols = dim(theta_samples)[2] - 1 
  
  eta_samples_matrix = matrix(nrow = samp_size, ncol = num_etas)
  
  #For each mcmc run
  for(i in 1:samp_size){
    R0 = theta_samples[i, 1]; k = mcmc_samples[i, 2]
    
    #print(paste0('dim2: ', length(theta_samples[i, 3:dim_cols])))
    density_samples = dgamma(theta_samples[i, 3:dim_cols], 
                             shape = epidemic_data[1:num_etas]*k,
                             scale = R0*k, log = TRUE)
    
    eta_samples_matrix[i, ] = density_samples #dgamma(theta_samples[i,3:dim_cols], shape = epidemic_data*k, scale = R0X*k, log = TRUE)
  }
  
  return(eta_samples_matrix)
}


#********************************************************************
#*
#1. GET IMPORTANCE SAMPLING PROPOSAL (LOG)
#*
#********************************************************************

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

#************************
#MULTI_DIMENSIONAL PROPOSAL
#************************
GET_LOG_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                                         n_samples, dof = 3, prob = 0.95, 
                                         priors_sseb = list(exp_prior = c(1,0)),
                                         priors_ssnb = list(pk_ga_shape = 0.001, pk_ga_rte = 0.001, pr0_unif = c(1.0,4),
                                                            p_prob_unif = c(0,1), ),
                                      priors_ssir = list(pk_exp = c(1,0), pR0_exp = c(1,0))){               
  
  #PARAMETERS REQUIRED
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; 
  samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR
  means = colMeans(mcmc_samples)
  print(paste0('dim of means = ', dim(means)))
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal) 
  
  #PRIORS
  if(FLAGS_MODELS$SSEB){
    n_dim = dim(mcmc_samples)[2] #3 #PUT IN DICTIONARY
    theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = n_dim) 
    
  } else if (FLAGS_MODELS$SSNB){
    n_dim =  dim(mcmc_samples)[2] #2 
    theta_samples_prior = matrix(c(rgamma(samp_size_prior, shape = priors_ssnb$pk_ga_shape, rate = priors_ssnb$pk_ga_rte),
                                   runif(samp_size_prior,  min = priors_ssnb$pr0_unif[1], max = priors_ssnb$pr0_unif[2])),
                                 ncol = n_dim)
  } else if (FLAGS_MODELS$SSIR) {
    n_dim = dim(mcmc_samples)[2]
    param_priors = cbind(rexp(samp_size_prior, rate = priors_ssir$pk_exp[1]),
                                   rexp(samp_size_prior, rate = priors_ssir$pR0_exp[1]))
    
    eta_priors_matrix = GET_SAMPLES_ETA_PRIORS(param_priors, epidemic_data, samp_size_prior)
    
    theta_samples_prior = matrix(c(param_priors, eta_priors_matrix), ncol = n_dim)
    
  }
  
  #Proposal samples
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(means, each = n_samples), ncol = n_dim)
  
  log_proposal = dmvt(theta_samples - matrix_means,
                      sigma = cov(mcmc_samples), df = dof) #log = TRUE log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #PRIOR DENSITIES 
  if(FLAGS_MODELS$SSEB){
    log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE) 
  } else if (FLAGS_MODELS$SSNB){
    log_priors = dgamma(theta_samples[,1], shape = priors_ssnb$pk_ga_shape, rate = priors_ssnb$pk_ga_rte, log = TRUE) +
      dunif(theta_samples[, 2], min = priors_ssnb$pr0_unif[1], max = priors_ssnb$pr0_unif[2], log = TRUE)
  } else if (FLAGS_MODELS$SSIR){

    log_density_eta_priors = GET_DENSITY_ETA_PRIORS(theta_samples, epidemic_data)
    log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) + log_density_eta_priors
  }
  
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_priors)) #1 x n_samples
  
  #max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_priors)
  #log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_priors - max_el))
  #log_q_s = LOG_SUM_EXP(log_q2) #LOG SUM EXP OF TWO COMPONENTS - See if the same
  #print(paste0('log_q_s2 ', log_q_s))
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_priors = log_priors)
  
  return(imp_samp_comps)  
}


#**************************************************************************************
#*
#* 2. GET P_HATS ESTIMATE OF MODEL EVIDENCE (LOG)
#*
#*****************************************************************************************
GET_LOG_MODEL_EVIDENCE_BASELINE <- function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
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
#* 2b. ESTIMATE OF MODEL EVIDENCE (PHAT) FOR SSEB MODEL  
#*******************************************************************
GET_LOG_MODEL_EVIDENCE <- function(mcmc_samples, epidemic_data,
                          FLAGS_MODELS, n_samples = 1000,
                          priors_ssnb = list(pk_ga_shape = 0.001, pk_ga_rte = 0.001, pr0_unif = c(1.0,4),
                                             p_prob_unif = c(0,1))){    
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(NA, n_samples)
  lambda_vec = get_lambda(epidemic_data); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_priors = imp_samp_comps$log_priors
  
  #SSEB MODEL 
  if (FLAGS_MODELS$SSEB){
    
    #GET ESTIMATE
    for (i in 1:n_samples) {
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3])
      
      vector_estimate_terms[i] = loglike + log_priors[i] - log_q[i]
    }
    
    #SSNB MODEL
  } else if (FLAGS_MODELS$SSNB) {
    
    for (i in 1:n_samples) {
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_SSNB(epidemic_data, lambda_vec, theta_samples[i,]) 
      
      vector_estimate_terms[i] = loglike + log_priors[i] - log_q[i]
    }
  } else if (FLAGS_MODELS$SSIR) {
    
    infectivity_vec = get_infectious_curve(epidemic_data)
    num_etas = length(epidemic_data)-1
    
    for (i in 1:n_samples) {
      loglike = LOG_LIKE_SSIR(epidemic_data, infectivity_vec, theta_samples[i, 1:2],  theta_samples[i, 3:(2+num_etas)]) 
      
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
LOAD_MCMC_GET_MODEL_EVIDENCE <- function(epidemic_data, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                start = 1, burn_in_pc = 0.2, BURN_IN = TRUE,
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
      
      print(paste0('i = ', i))
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
      
      print(paste0('i = ', i))
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
      
      print(paste0('i = ', i))
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
    
  }
  
  return(estimates_vec) 
}
