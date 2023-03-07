#Model Comparison - Importance Sampled
#'Model evidence estimator via importance sampling'

#LIBRARIES
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
#library(compositions)

#***************************************
#*
#1. GET IMPORTANCE SAMPLING PROPOSAL
#*
#***************************************
GET_LOG_Q_PROPOSAL_UNI_VAR <- function(mcmc_samples, epidemic_data, #priors = 
                                 n_samples, num_dims = 1) {               #1.OTHER: GET SINGLE DIM PROPSAL. SINGLE T DISTRIBUTION
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  samp_size_proposal = round(0.95*n_samples); samp_size_prior =  n_samples - samp_size_proposal
  
  #*******
  #THETA SAMPLES: PROPOSAL + PRIOR
  mean_mcmc = mean(mcmc_samples)
  theta_samples_proposal = rt(samp_size_proposal, df = num_dims) + mean_mcmc 
  theta_samples_prior = c(rexp(samp_size_prior))
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  proposal = dt(theta_samples - mean_mcmc, df = num_dims, log = FALSE)
  prior = dexp(theta_samples[,1])
  q = 0.95*proposal + 0.05*prior
  log_q = LOG_SUM_EXP(q) #LOG SUM EXP OF TWO COMPONENTS
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps)
}

#MUTLI DIM PROPOSAL
GET_LOG_Q_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data,  #GET_PROPOSAL_MULTI_DIM
                                   n_samples) {               
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  samp_size_proposal = round(0.95*n_samples); samp_size_prior =  n_samples - samp_size_proposal
  
  #THETA SAMPLES: PROPOSAL + PRIOR
  means = colMeans(mcmc_samples)
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = 3) + means 
  theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), rexp(samp_size_prior)), ncol = 3) 
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  proposal = dmvt(theta_samples - means, sigma = cov(mcmc_samples), df = 3, log = FALSE)
  prior = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1))
  q = 0.95*proposal + 0.05*prior
  log_q = LOG_SUM_EXP(q) #LOG SUM EXP OF TWO COMPONENTS
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q)
  
  return(imp_samp_comps) #log_q 
}

#*********************************************************
#*
#* 2. GET LOG P_HATS 
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


#MULTI PARAM MODELS 
GET_LOG_P_HAT <-function(mcmc_samples, epidemic_data, 
                                     FLAGS_LIST = list(SSEB = TRUE,
                                                       SSIB = FALSE, SSIC = FALSE),
                                     n_samples = 10000) {
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  sum_estimate = 0
  imp_samp_comps = GET_LOG_Q_PROPOSAL_MULTI_DIM(mcmc_samples, epidemic_data, n_samples)
  theta_samples = imp_samp_comps$theta_samples #Some samples could be negative
  log_q = imp_samp_comps$log_q
  
  #PRIORS 
  if (FLAGS_LIST$SSEB | FLAGS_LIST$SSIB) {
    priors = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) + dexp((theta_samples[,3] - 1))
    lambda_vec = get_lambda(epidemic_data); 
    
  } else {
    priors = dexp(theta_samples[,1]) + dexp(theta_samples[,2]) 
    infectivity = get_infectious_curve(epidemic_data)
  }
  
  #LOG SUM EXP (LOOP)
  vector_log_sum_exp = c()
  for(i in 1:n_samples){
    
    if(FLAGS_LIST$SSEB) {
      loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],  theta_samples[i, 2],
                              theta_samples[i, 3])
      if (is.na(loglike)){
        vector_log_sum_exp[i] = log(priors[i]) - log_q
      } else {
        vector_log_sum_exp[i] = loglike + log(priors[i]) - log_q
      }
      
      vector_log_sum_exp[i] =  + log(priors[i]) - log_q
      
    } else if(FLAGS_LIST$SSIB) {
      
      loglike = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
                             theta_samples[i, 3])
      if (is.na(loglike)){
        vector_log_sum_exp[i] = log(priors[i]) - log_q
      } else {
        vector_log_sum_exp[i] = loglike + log(priors[i]) - log_q
      }
      
      vector_log_sum_exp[i] =  + log(priors[i]) - log_q
      
    } else if (FLAGS_LIST$SSIC) {
      
      vector_log_sum_exp[i] = LOG_LIKE_SSI(epidemic_data, theta_samples[i, 1],  theta_samples[i, 2],
                              theta_samples[i, 3]) + log(priors[i]) - log_q
    }
  }

  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_log_sum_exp)
  
  return(log_p_hat)
}

#SSEB
phat_sseb = GET_LOG_P_HAT(mcmc_samples, data_baseI)
#Base
phat_base = GET_LOG_P_HAT_BASELINE(mcmc_samples, data_baseI)
#SSIB
phat_ssib = GET_LOG_P_HAT(mcmc_samples, data_baseI, FLAGS_LIST = list(SSEB = FALSE, SSIB = TRUE,
                                                                      SSIC = FALSE))

#*********************************************************
#*
#* 3. GET POSTERIOR MODEL PROBABILITIES
#* 
#* **********************************************************
GET_POSTERIOR_MODEL_PROB <- function(num_models = 3, 
                                     probs_models = list(prob_mech1 = 1/3, prob_mech2 = 1/3), #mech = mechanism; baseline (0.5) or sse (0.25)
                                     log_phats = list(mod1 = mod1,
                                                           mod2 = mod2, mod3 = mod3)){ 
  
  #PARAMS
  vec_model_diffs = c()
  
  for(i in 2:num_models){
    print(paste0('i = ', i))
    vec_model_diffs[i-1] = exp(log(probs_models[[2]]) + log_phats[[i]] -
                                 log(probs_models[[1]]) - log_phats[[1]])
  }
  
  posterior_prob =  1/(1 + sum(vec_model_diffs))
  
  return(posterior_prob)
}

# GET_AGGREGATE_POSTERIOR_MODEL_PROB
GET_AGG_POSTERIOR_MODEL_PROB <- function(num_models = 3, 
                                               probs_models = list(prob_mech1 = 1/3, prob_mech2 = 1/3), #mech = mechanism; baseline (0.5) or sse (0.25)
                                               list_log_phats = list(mod1 = mod1,
                                                                     mod2 = mod2, mod3 = mod3)){ 
  
  #FOR EACH REP
  num_reps = length(list_log_phats$mod1)
  vec_post_probs = c()
  
  for (i in 1:num_reps){
  
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(num_models = 3, 
                                       probs_models = probs_models,
                                       log_phats = list(mod1 = list_log_phats$mod1[i],
                                                        mod2 = list_log_phats$mod2[i],
                                                        mod3 = list_log_phats$mod3[i]))
    
  }

  return(vec_post_probs)
  
}

#APPLY
post_prob1 = GET_POSTERIOR_MODEL_PROB(log_phats = list(mod1 = phat_base,
                                                       mod2 = phat_sseb, mod3 = phat_ssib))

post_prob2 = GET_POSTERIOR_MODEL_PROB(log_phats = list(mod1 = phat_sseb,
                                                       mod2 = phat_base, mod3 = phat_ssib))
post_prob2

#APPLY AGGREGATE
vec_post_probs = GET_AGG_POSTERIOR_MODEL_PROB(list_log_phats = list(mod1 = c(-101, -102, -103, 102, 101),
                                                         mod2 = c(-113, -114, 112, 111, -115), mod3 = c(-15, -15.5, -16, -16.1, -15.9) ))
vec_post_probs

#******************************************************************************
#*
# 4. LOAD MCMC + GET POSTERIOR MODEL PROBABILITIES
#*
#******************************************************************************

LOAD_MCMC_GET_P_HAT <- function(epidemic_data, OUTPUT_FOLDER, run = 1, n_repeats = 100,
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
      log_phat = GET_LOG_P_HAT_BASELINE(mcmc_samples$r0_vec, epidemic_data) 
      estimates_vec[i] = log_phat
      print(estimates_vec)
    }
    #SAVE ESTIMATES
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_base_vec', i ))
    
  } else if(FLAGS_LIST$SSEB){
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
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_sseb_vec', i ))
    
  } else if (FLAGS_LIST$SSIB){
    
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
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_ssib_vec', i ))
    
  }
  
  return(estimates_vec) 
}



#RUN MCMC + GET ESTS
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
