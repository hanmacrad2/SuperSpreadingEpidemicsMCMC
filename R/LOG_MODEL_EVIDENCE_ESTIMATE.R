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
GET_LOG_Q_PROPOSAL_UNI_VAR <- function(mcmc_samples, epidemic_data, 
                                 n_samples, num_dims = 1) {               
  
  #Code
  'Fix samp_proposal'
  
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
  theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = 3) 
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DEFENSE MIXTURE
  proposal = dmvt(theta_samples - means, sigma = cov(mcmc_samples), df = 3, log = FALSE)
  prior = dexp(theta_samples[,1])*dexp(theta_samples[,2])*dexp((theta_samples[,3] - 1))
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


#MULTI PARAMETER MODELS 
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
#***********************************************************
GET_POSTERIOR_MODEL_PROB <- function(num_models = 3, 
                                     probs_models = list(prob_mech1 = 0.25, prob_mech2 = 0.5, prob_mech3 = 0.25), #mech = mechanism; baseline (0.5) or sse (0.25)
                                     log_phats = list(mod1 = mod1,
                                                      mod2 = mod2, mod3 = mod3)){ 
  
  #PARAMS
  vec_model_diffs = c()
  
  for(i in 1:num_models-1){
    vec_model_diffs[i] = exp(log(probs_models[[i+1]]) + log_phats[[i+1]] -
                                 log(probs_models[[1]]) - log_phats[[1]])
  }
  
  print(paste0('sum(vec_model_diffs) = ', sum(vec_model_diffs)))
  
  posterior_prob =  1/(1 + sum(vec_model_diffs))
  
  print(paste0('posterior_prob= ', posterior_prob))
  
  return(posterior_prob)
}

#GET_AGGREGATE_POSTERIOR_MODEL_PROB
GET_AGG_POSTERIOR_PROB <- function(num_models = 3, 
                                   probs_models = list(prob_mech1 = 0.25, prob_mech2 = 0.5, prob_mech3 = 0.25),
                                               list_log_phats = list(mod1 = mod1,
                                                                     mod2 = mod2, mod3 = mod3)){ 
  'Get posterior probabilites for mulitple P_hat reps'
  
  #FOR EACH REP
  num_reps = length(list_log_phats$mod1)
  vec_post_probs = c()
  
  for (i in 1:num_reps){
    print(paste0('rep = ', i))
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(num_models = 3, 
                                       probs_models = probs_models,
                                       log_phats = list(mod1 = list_log_phats$mod1[i],
                                                        mod2 = list_log_phats$mod2[i],
                                                        mod3 = list_log_phats$mod3[i]))
    
  }

  return(vec_post_probs)
  
}

#*********************
#* 4. GET BAYES FACTORS
#*********************
GET_LOG_BAYES_FACTORS <- function (list_log_phat_mod1, list_log_phat_mod2){
  
  log_bf = list_log_phat_mod1 - list_log_phat_mod2
  
  return(log_bf)
  
}

#******************************************************************************
#*
# 5. LOAD MCMC + GET POSTERIOR MODEL PROBABILITIES
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
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_base_', run, '.rds' ))
    
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
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_sseb_', run, '.rds' ))
    
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
    saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_ssib_', run, '.rds' ))
    
  }
  
  return(estimates_vec) 
}

#******************
#**NEED TO EDIT
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

#******************************************************************************
#* PLOTTING RESULTS FUNCTION
#*
#******************************************************************************
PLOT_MODEL_COMPARISON_RESULTS <- function(model_comp_results, 
                                          result_type = 'Bayes Factors: Baseline vs SSEB Models. ',
                                          data_type = 'Baseline', 
                                  n_reps = 100, FLAG_RESULT_TYPE = list(log = TRUE)){
  
  #TITLE
  if (FLAG_RESULT_TYPE$log) {
    #axis_label = paste0(result_type, ' (log).')
    axis_label = paste0(result_type)
  } else  axis_label = paste0(result_type)
  
  #Title
  titleX = paste0(axis_label, data_type, ' data. ', n_reps, ' reps.')
  
  #PLOT
  par(mfrow = c(2,1))
  boxplot(model_comp_results,
          ylab = axis_label,
          main = titleX)
  
  hist(model_comp_results, breaks = 100, freq = FALSE,
       xlab = axis_label,
       main = titleX)
  
}


#MODEL EVIDENCE RESULTS
PLOT_MODEL_EV_RESULTS <- function(posterior_results, model_type = 'Baseline', data_type = 'Baseline', 
                                  n_reps = 100, FLAG_RESULT_TYPE = list(phat = FALSE, post_prob = TRUE,
                                                                        log = TRUE)){
  
  #TITLE
  if(FLAG_RESULT_TYPE$phat) result_type = 'P hat, '
  if(FLAG_RESULT_TYPE$post_prob) result_type = 'Posterior model probability '
  #LOG = TRUE
  if (FLAG_RESULT_TYPE$log) {
    axis_label = paste0(result_type, '(log), ', model_type, ' model. ')
    posterior_results = log(posterior_results)
  } else  axis_label = paste0(result_type, model_type, ' model. ')
  
  #Title
  titleX = paste0(axis_label, data_type, ' data. ', n_reps, ' reps.')

  #PLOT
  par(mfrow = c(2,1))
  boxplot(posterior_results,
          ylab = axis_label,
          main = titleX)
  
  hist(posterior_results, breaks = 50, freq = FALSE,
       xlab = axis_label,
       main = titleX)
  
}

#*******************
#* APPLICATION OF FUNCTIONS
#APPLY
post_prob1 = GET_POSTERIOR_MODEL_PROB(log_phats = list(mod1 = phat_base,
                                                       mod2 = phat_sseb, mod3 = phat_ssib))

post_prob2 = GET_POSTERIOR_MODEL_PROB(log_phats = list(mod1 = phat_sseb,
                                                       mod2 = phat_base, mod3 = phat_ssib))

#APPLY AGGREGATE
vec_post_probs = GET_AGG_POSTERIOR_MODEL_PROB(list_log_phats = list(mod1 = c(-101, -102, -103, 102, 101),
                                                                    mod2 = c(-113, -114, 112, 111, -115), mod3 = c(-15, -15.5, -16, -16.1, -15.9) ))

