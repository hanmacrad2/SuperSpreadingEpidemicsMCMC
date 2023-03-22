#'Model evidence estimator via importance sampling'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
#library(compositions)

#***************************************
#*
#1. GET IMPORTANCE SAMPLING PROPOSAL (LOG)
#*
#***************************************
GET_LOG_Q_PROPOSAL_UNI_VAR <- function(mcmc_samples, epidemic_data, 
                                       n_samples, dof = 3, prob = 0.95) {   
  
  'Get proposal q for univariate dist'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  sampling_prob = rbinom(n_samples, 1, prob)
  samp_size_proposal = length(which(sampling_prob == 1)); samp_size_prior =  n_samples - samp_size_proposal
  print(paste0('samp_size_proposal = ', samp_size_proposal)); print(paste0('samp_size_prior = ', samp_size_prior))
  prob_prop = samp_size_proposal/n_samples; prob_prior = 1 - prob_prop
  
  #*******
  #THETA SAMPLES: PROPOSAL + PRIOR
  mean_mcmc = mean(mcmc_samples); sd_mcmc = sd(mcmc_samples)
  theta_samples_proposal = sd_mcmc*rt(samp_size_proposal, df = dof) + mean_mcmc 
  theta_samples_prior = c(rexp(samp_size_prior))
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dt(theta_samples - mean_mcmc, df = dof, log = TRUE) - log(sd_mcmc)
  log_prior = dexp(theta_samples[,1], log = TRUE)
  
  #LOG SUM EXP TRICK TO GET LOG_Q
  max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_prior)
  log_q = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_prior - max_el))
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps)
}

#MUTLI DIM PROPOSAL
GET_LOG_Q_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data,  #GET_PROPOSAL_MULTI_DIM
                                         n_samples, dof = 3, prob = 0.95) {               
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  # sampling_prob = rbinom(n_samples, 1, prob)
  # samp_size_proposal = length(which(sampling_prob == 1)); samp_size_prior =  n_samples - samp_size_proposal
  # print(paste0('samp_size_proposal = ', samp_size_proposal)); print(paste0('samp_size_prior = ', samp_size_prior))
  # prob_prop = samp_size_proposal/n_samples; prob_prior = 1 - prob_prop
  
  #SAMPLING SIZE #2
  samp_size_proposal = prob*n_samples; samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1-prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR
  means = colMeans(mcmc_samples)
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) + means 
  theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = 3) 
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dmvt(theta_samples - means, sigma = cov(mcmc_samples), df = dof, log = TRUE)
  log_prior = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) + dexp((theta_samples[,3] - 1), log = TRUE)
  
  max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_prior)
  log_q = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_prior - max_el))
  
  #log_q_s = LOG_SUM_EXP(log_q) #LOG SUM EXP OF TWO COMPONENTS
  
  #print(paste0('LOG_Q: ', log_q))
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps)  
}

#*********************************************************
#*
#* 2. GET P_HATS ESTIMATE OF MODEL EVIDENCE (LOG)
#*
#************************************************************
GET_LOG_P_HAT_BASELINE <-function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0
  imp_samp_comps = GET_LOG_Q_PROPOSAL_UNI_VAR(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples 
  log_q = imp_samp_comps$log_q
  
  #PRIORS 
  priors = dexp(theta_samples) 
  
  #LOG SUM EXP (LOOP)
  vector_log_sum_exp = c()
  for(i in 1:n_samples){         
    
    loglike = LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])
    
    #NEGATIVE THETA SAMPLES -> NEGATIVE LOGLIKELIHOOD
    if (is.na(loglike)){
      vector_log_sum_exp[i] = log(priors[i]) - log_q
    } else {
      vector_log_sum_exp[i] = loglike + log(priors[i]) - log_q
    }
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
  
  return(log_p_hat)
}


#MULTI PARAMETER MODELS 
GET_LOG_P_HAT <- function(mcmc_samples, epidemic_data, 
                          FLAGS_MODELS = list(SSEB = TRUE,
                                              SSIB = FALSE, SSIC = FALSE),
                          n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0
  imp_samp_comps = GET_LOG_Q_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples 
  log_q = imp_samp_comps$log_q
  
  #PRIORS 
  if (FLAGS_MODELS$SSEB | FLAGS_MODELS$SSIB) {
    print(paste0('FLAGS_MODELS$SSEB' = FLAGS_MODELS$SSEB))
    log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE)
    lambda_vec = get_lambda(epidemic_data); 
    
  } else {
    log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) 
    infectivity = get_infectious_curve(epidemic_data)
  }
  
  #LOG SUM EXP (LOOP)
  vector_log_sum_exp = c()
  for(i in 1:n_samples){
    
    if(i%%100 == 0) print(i)
    
    if(FLAGS_MODELS$SSEB) {
      loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
                              theta_samples[i, 3])
      if (is.na(loglike)){
        vector_log_sum_exp[i] = log_priors[i] - log_q
      } else {
        vector_log_sum_exp[i] = loglike + log_priors[i] - log_q
        #print(paste0('vector_log_sum_exp[i]', vector_log_sum_exp[i]))
      }
      
    } #else if(FLAGS_MODELS$SSIB) {
    
    # loglike = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
    #                        theta_samples[i, 3])
    # if (is.na(loglike)){
    #   vector_log_sum_exp[i] = log(priors[i]) - log_q
    # } else {
    #   vector_log_sum_exp[i] = loglike + log(priors[i]) - log_q
    # }
    # 
    # vector_log_sum_exp[i] =  + log(priors[i]) - log_q
    
    # } else if (FLAGS_MODELS$SSNB) {
    #   
    #   vector_log_sum_exp[i] = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
    #                           #theta_samples[i, 3]) + log(priors[i]) - log_q
    # }
  }
  
  #print(paste0('LOG_SUM_EXP(vector_log_sum_exp)',  LOG_SUM_EXP(vector_log_sum_exp)))
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
  
  return(log_p_hat)
}

#******************************************************************************
#*
# 3. LOAD MCMC + GET P_HAT ESTIMATES MODEL EVIDENCE ESTIMATES
#*
#******************************************************************************
LOAD_MCMC_GET_P_HAT <- function(epidemic_data, OUTPUT_FOLDER, run = 1, n_repeats = 100,
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                    SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #FOLDER
  CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)
  create_folder(CURRENT_OUTPUT_FOLDER)
  
  #Parameters
  estimates_vec = c()
  
  if (FLAGS_MODELS$BASE){
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #READ SAMPLES
      mcmc_samples = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base_', i ))
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      log_phat = GET_LOG_P_HAT_BASELINE(mcmc_samples$r0_vec, epidemic_data) 
      estimates_vec[i] = log_phat
      print(estimates_vec)
    }
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_base_', run, '.rds' ))
    
  } else if(FLAGS_MODELS$SSEB){
    print('sseb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      mcmc_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', i ))
      mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data) 
      
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_sseb_', run, '.rds' ))
    
  } else if (FLAGS_MODELS$SSIB){
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      mcmc_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssib_', i ))
      mcmc_samples =  matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data) 
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_ssib_', run, '.rds' ))
    
  }
  
  return(estimates_vec) 
}
