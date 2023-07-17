#Model Comparison - Importance Sampled
#'Model evidence estimator via importance sampling'

#LIBRARIES
library(SuperSpreadingEpidemicsMCMC)
library(compositions)
#library(mvtnorm)

#***************************************
#1. GET IMPORTANCE SAMPLING PROPOSAL
#***************************************

#*****************
GET_PROPOSAL_UNI_VAR <- function(mcmc_samples, epidemic_data, #priors = 
                                 n_samples) {               #1.OTHER: GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #PROPOSAL 
  means = mean(mcmc_samples)
  theta_samples = rlnorm(n_samples, log(mean(mcmc_samples)), sd(mcmc_samples))
  proposal = dlnorm(theta_samples, log(mean(mcmc_samples)), sd(mcmc_samples))
  imp_samp_comps = list(theta_samples = theta_samples, proposal = proposal)
  
  return(imp_samp_comps)
}

GET_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data, #priors = 
                                   n_samples) {               #1.OTHER: GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #PROPOSAL 
  means = colMeans(mcmc_samples)
  theta_samples = rlnorm.rplus(n_samples, log(colMeans(mcmc_samples)), cov(mcmc_samples))
  proposal = dlnorm.rplus(theta_samples, log(colMeans(mcmc_samples)), cov(mcmc_samples))
  imp_samp_comps = list(theta_samples = theta_samples, proposal = proposal)
  
  #t dist (incompatible as not strictly positive)
  #theta_samples = rmvt(n_samples, sigma = cov(mcmc_samples), df = 3) + means #Samples (theta)
  #proposal = dmvt(theta_samples - means, sigma = cov(mcmc_output), df = 3, log = FALSE)
  
  return(imp_samp_comps)
}

#***************************************************
#* GET MODEL EVIDENCE ESTIMATES VIA IMPORTANCE SAMPLING
#*****************************************************

#**********************************************
#* 1. PHAT BASE MODEL
GET_IMP_SAMP_MODEL_EV_BASE <-function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0
  imp_samp_comps = GET_PROPOSAL_UNI_VAR(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  proposal = imp_samp_comps$proposal
  
  #PRIORS 
  priors = dexp(theta_samples)  #dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    
    estimate = LOG_LIKE_BASELINE(epidemic_data, theta_samples[i]) +
      log(priors[i]) - log(q_defense_mixture[i])
    
    if(is.infinite(estimate)){
      print('yes infinite')
      print(paste0('LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])', LOG_LIKE_BASELINE(epidemic_data, theta_samples[i])))
      print(paste0('log(priors[i])', log(priors[i])))
      print(paste0('theta_samples[i]', theta_samples[i]))
      print(paste0('log(q_defense_mixture[i])', log(q_defense_mixture[i])))
    }
    
    #print(estimate)
    sum_estimate = sum_estimate + estimate
  }
  
  p_hat_est = sum_estimate/n_samples
  
  return(p_hat_est)
}

#**********************************************
#* 2. PHAT SSEB 
GET_IMP_SAMP_MODEL_EV_SSEB <-function(mcmc_samples, epidemic_data, n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data); sum_estimate = 0
  imp_samp_comps = GET_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  proposal = imp_samp_comps$proposal
  
  #PRIORS 
  priors = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1)) #LOGS
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    
    estimate = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
                             theta_samples[i, 3]) +
      log(priors[i]) - log(q_defense_mixture[i])
    
    #print(estimate)
    sum_estimate = sum_estimate + estimate
  }
  
  p_hat_est = sum_estimate/n_samples
  
  return(p_hat_est)
}

#***************
#* PHAT SSIB

GET_IMP_SAMP_MODEL_EV_SSB <-function(mcmc_samples, epidemic_data, 
                                     FLAGS_LIST = list(SSEB = FALSE,
                                                       SSIB = TRUE, SSIC = FALSE),
                                     n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0
  imp_samp_comps = GET_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  proposal = imp_samp_comps$proposal
  
  #PRIORS 
  if (FLAGS_LIST$SSEB | FLAGS_LIST$SSIB) {
    priors = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1))
    lambda_vec = get_lambda(epidemic_data); 
    
  } else {
    priors = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) 
    infectivity = get_infectious_curve(epidemic_data)
  }
  
  #MIXTURE
  q_defense_mixture = 0.95*proposal + 0.05*priors
  
  #LOOP
  for(i in 1:n_samples){
    
    if(FLAGS_LIST$SSEB) {
      
      estimate = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
                               theta_samples[i, 3])
      
    } else if(FLAGS_LIST$SSIB) {
      
      estimate = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
                              theta_samples[i, 3])
      
    } else if (FLAGS_LIST$SSIC) {
      
      estimate = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
                              theta_samples[i, 3])
    }
    
    sum_estimate = sum_estimate + estimate + log(priors[i]) - log(q_defense_mixture[i])
  }
  
  p_hat_est = sum_estimate/n_samples
  
  return(p_hat_est)
}


#******************************************************************************
# 1. RUN AUTOMATICALLY 
#******************************************************************************
RUN_MCMC_MODEL_EV_IMP_SAMP <- function(epidemic_data, OUTPUT_FOLDER, run = 1, n_repeats = 100,
                                       FLAGS_LIST = list(BASE = FALSE, SSEB = FALSE,
                                                         SSIB = TRUE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #FOLDER
  CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)
  create_folder(CURRENT_OUTPUT_FOLDER)
  
  #Parameters
  estimates_vec = c()
  
  if (FLAGS_LIST$BASE){
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #MCMC SAMPLES
      mcmc_samples = MCMC_INFER_BASELINE(epidemic_data)
      #SAVE MCMC
      saveRDS(mcmc_samples, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base_', i ))
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_IMP_SAMP_MODEL_EV_BASE(mcmc_samples$r0_vec, epidemic_data) 
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
    }
  } else {
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      
      #MCMC SAMPLES
      if(FLAGS_LIST$SSEB){
        mcmc_output = MCMC_INFER_SSEB(epidemic_data)
        saveRDS(mcmc_output, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', i ))
        mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)
        
      } else if (FLAGS_LIST$SSIB){
        mcmc_output = MCMC_INFER_SSIB(epidemic_data)
        saveRDS(mcmc_output, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssib_', i ))
        mcmc_samples =  matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
        
      } else if (FLAGS_LIST$SSIC){
        mcmc_output = MCMC_INFER_SSIC(epidemic_data)
        saveRDS(mcmc_output, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssic_', i )) 
        #SSIC_PARAMS + ETA
        #mcmc_samples = mcmc_output$
        
      } else if (FLAGS_LIST$SSEC) {
        mcmc_output = SSI_MCMC_ADAPTIVE(epidemic_data)
        saveRDS(mcmc_output, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssec_', i )) 
      }
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_IMP_SAMP_MODEL_EV_SSB(mcmc_samples, epidemic_data) 
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
    }
    
  }
  return(estimates_vec) 
}

#FOLDER SAVE
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/BASE"

#APPLY
RUN_MCMC_MODEL_EV_IMP_SAMP(epidemic_data, OUTPUT_FOLDER)


#**************************
#* LOAD_MCMC_GET_P_EST

LOAD_MCMC_GET_P_EST <- function(epidemic_data, OUTPUT_FOLDER, run = 1, n_repeats = 100,
                                FLAGS_LIST = list(BASE = FALSE, SSEB = TRUE,
                                                  SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #FOLDER
  CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)
  create_folder(CURRENT_OUTPUT_FOLDER)
  
  #Parameters
  estimates_vec = c()
  
  if (FLAGS_LIST$BASE){
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #READ SAMPLES
      mcmc_samples = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base_', i ))
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_IMP_SAMP_MODEL_EV_BASE(mcmc_samples$r0_vec, epidemic_data) 
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
    }
    #SAVE ESTIMATES
    
  } else if(FLAGS_LIST$SSEB){
    print('sseb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      mcmc_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', i ))
      mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_IMP_SAMP_MODEL_EV_SSB(mcmc_samples, epidemic_data) 
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_sseb_vec', i ))
    
  } else if (FLAGS_LIST$SSIB){
    
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      mcmc_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssib_', i ))
      mcmc_samples =  matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
      
      #GET PHAT ESTIMATE OF MODEL EVIDENCE
      phat_estimate = GET_IMP_SAMP_MODEL_EV_SSB(mcmc_samples, epidemic_data) 
      estimates_vec[i] = phat_estimate
      print(estimates_vec)
      
    }
    
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_ssib_vec', i ))
    
  }
  
  return(estimates_vec) 
}