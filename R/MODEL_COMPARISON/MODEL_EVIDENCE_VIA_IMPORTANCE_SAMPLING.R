#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
#library(compositions)

#*************************************
#* LOG SUM EXP
#* ***********************************
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
GET_LOG_PROPOSAL_Q_UNI_VAR <- function(mcmc_samples, epidemic_data, 
                                       n_samples, dof = 3, prob = 0.95) {   
  
  'Get proposal q for univariate dist'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  
  #SAMPLING SIZE
  samp_size_proposal = prob*n_samples; samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1-prob_prop
  
  #*******
  #THETA SAMPLES: PROPOSAL + PRIOR
  mean_mcmc = mean(mcmc_samples); sd_mcmc = sd(mcmc_samples)
  theta_samples_proposal = sd_mcmc*rt(samp_size_proposal, df = dof) + mean_mcmc 
  theta_samples_prior = c(rexp(samp_size_prior))
  #theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  theta_samples = c(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  log_proposal = dt((theta_samples - mean_mcmc)/sd_mcmc, df = 1, log = TRUE) - log(sd_mcmc) #ADDED
  
  print(paste0('1. mean log_proposal = ', mean(log_proposal)))
  
  #dmvt: not matching ***
  log_proposal2 = dmvt(matrix(theta_samples),
                       sigma = diag(sd_mcmc^2, 1), df = 1, log = TRUE)
  
  log_proposal2 = dmvt(matrix((theta_samples - mean_mcmc)/sd_mcmc,
                               rep(1, length(theta_samples)), rep(1, length(theta_samples))), 
                       sigma = diag(sd_mcmc, length(theta_samples)), log=TRUE)
  print(paste0('2. mean log_proposal2 = ', mean(log_proposal2)))
  #print('')
  
  #dmvt(matrix(1,1,1), sigma=matrix(4,1,1), log=TRUE)
  #dt((1-0)/2, df=1, log=TRUE)-log(2)
  
  log_prior = dexp(theta_samples, log = TRUE) #CHECK [,1]
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_prior)) #CALCULATE WITH LOG SUM EXP TRICK ASWELL & SEE IF MATCH
  print(paste0('mean log_q = ', mean(log_q)))
  
  #LOG SUM EXP TRICK TO GET LOG_Q (MATCH) 
  max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_prior)
  log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_prior - max_el))
  print(paste0('4 mean log_q2 lse = ', mean(log_q2)))
  #print('')
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps)
}

#*************************************
#1b. MUTLI DIM PROPOSAL Q
#*************************************
GET_LOG_PROPOSAL_Q_MULTI_DIM <- function(mcmc_samples, epidemic_data,  
                                         n_samples, dof = 3, prob = 0.95) {               
  
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
  
  theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = 3) 
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  print(paste0('mean_mcmc = ', means))
  #print(paste0('mean proposal = ', colMeans(theta_samples_proposal)))
  #print(paste0('mean prior = ', colMeans(theta_samples_prior)))
  
  #DEFENSE MIXTURE
  #log_proposal = dmvt(theta_samples - means, sigma = cov(mcmc_samples), df = dof, log = TRUE)
  log_proposal = dmvt(theta_samples - matrix( rep(means, each = n_samples), ncol = 3),
                      sigma = cov(mcmc_samples), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  = 2, f(x,y) = 0.2) for exmaples
  log_prior = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) + dexp((theta_samples[,3] - 1), log = TRUE)
  
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_prior)) #1 x n_samples
  
  # print(paste0('mean_mcmc = ', means))
  # print(paste0('mean log_proposal = ', mean(log_proposal)))
  # print(paste0('mean prior = ', mean(log_prior, na.rm = TRUE)))
  # print(paste0('1 log_q = ', log_q))
  # 
  # max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_prior)
  # log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_prior - max_el))
  # log_q_s = LOG_SUM_EXP(log_q2) #LOG SUM EXP OF TWO COMPONENTS - See if the same
  # print(paste0('2 log_q_s = ', log_q_s))
  
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
  imp_samp_comps = GET_LOG_PROPOSAL_Q_UNI_VAR(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples 
  log_q = imp_samp_comps$log_q
  
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
#*******************************************************************
GET_LOG_P_HAT_SSEB <- function(mcmc_samples, epidemic_data, n_samples = 1000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0; lambda_vec = get_lambda(epidemic_data); 
  imp_samp_comps = GET_LOG_PROPOSAL_Q_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples ; log_q = imp_samp_comps$log_q
  
  #PRIORS 
  log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
    dexp((theta_samples[,3] - 1), log = TRUE)
  
  #LOG SUM EXP (LOOP)
  vector_terms = rep(NA, n_samples) 
  for(i in 1:n_samples){
    if(i%%100 == 0) print(i)
    
    loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
                            theta_samples[i, 3])
    vector_terms[i] = loglike + log_priors[i] - log_q[i]
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  log_p_hat2 = -log(n_samples) + log(sum(exp(vector_terms)))
  print(paste0('log_p_hat2 = ', log_p_hat2))
  
  return(log_p_hat)
}


#MULTI PARAMETER MODELS 
# GET_LOG_P_HAT <- function(mcmc_samples, epidemic_data, 
#                           FLAGS_MODELS = list(SSEB = TRUE,
#                                               SSIB = FALSE, SSIC = FALSE),
#                           n_samples = 1000) {
#   
#   'Estimate of model evidence for SSEB model using Importance Sampling'
#   
#   #PARAMS
#   sum_estimate = 0
#   imp_samp_comps = GET_LOG_PROPOSAL_Q_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
#   theta_samples = imp_samp_comps$theta_samples 
#   log_q = imp_samp_comps$log_q
#   
#   #PRIORS 
#   if (FLAGS_MODELS$SSEB | FLAGS_MODELS$SSIB) {
#     print(paste0('FLAGS_MODELS$SSEB' = FLAGS_MODELS$SSEB))
#     log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
#       dexp((theta_samples[,3] - 1), log = TRUE)
#     lambda_vec = get_lambda(epidemic_data); 
#     
#   } 
#   
#   #LOG SUM EXP (LOOP)
#   vector_log_sum_exp = rep(NA, n_samples) 
#   for(i in 1:n_samples){
#     
#     if(i%%100 == 0) print(i)
#     
#     if(FLAGS_MODELS$SSEB) {
#       loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
#                               theta_samples[i, 3])
#       if (is.na(loglike)){
#         vector_log_sum_exp[i] = log_priors[i] - log_q[i]
#       } else {
#         vector_log_sum_exp[i] = loglike + log_priors[i] - log_q[i]
#         #print(paste0('vector_log_sum_exp[i]', vector_log_sum_exp[i]))
#       }
#       
#     } #else if(FLAGS_MODELS$SSIB) {
#     
#     loglike = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
#                            theta_samples[i, 3])
#     if (is.na(loglike)){
#       vector_log_sum_exp[i] = log(priors[i]) - log_q
#     } else {
#       vector_log_sum_exp[i] = loglike + log(priors[i]) - log_q
#     }
# 
#     vector_log_sum_exp[i] =  + log(priors[i]) - log_q
# 
#     } else if (FLAGS_MODELS$SSNB) {
# 
#       vector_log_sum_exp[i] = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
#                               #theta_samples[i, 3]) + log(priors[i]) - log_q
#     }
#   }
#   
#   #print(paste0('LOG_SUM_EXP(vector_log_sum_exp)',  LOG_SUM_EXP(vector_log_sum_exp)))
#   log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
#   
#   return(log_p_hat)
# }


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
  estimates_vec = rep(NA, n_repeats) 
  
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
      #phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data) 
      phat_estimate = GET_LOG_P_HAT_SSEB(mcmc_samples, epidemic_data) 
      
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
