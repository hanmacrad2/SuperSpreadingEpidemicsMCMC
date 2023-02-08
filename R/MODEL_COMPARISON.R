#MODEL COMPARISON

#*********************
#* MODEL COMPARISON BY RATIO OF MODEL EVIDENCES
#***********

LOG_MODEL_EVIDENCE <- function(loglike_vec){
  
  'Model evidence via log-sum-exp trick'
  
  loglike_vec = - loglike_vec
  m = max(loglike_vec, na.rm = TRUE)
  log_model_ev = m + log(mean(exp(loglike_vec - m)))
  
  return(-log_model_ev)
}

#GET BAYES FACTORS
GET_BAYES_FACTORS <- function(loglike_vec1, loglike_vec2){
  
  'Get Bayes factor via ratio of the model evidence'
  
  #bayes_factor = MODEL_EVIDENCE(loglike_vec1)/MODEL_EVIDENCE(loglike_vec2) 
  log_bf = LOG_MODEL_EVIDENCE(loglike_vec1)  - LOG_MODEL_EVIDENCE(loglike_vec2)
  bayes_factor = exp(log_bf)
  
  return(bayes_factor)
}

#RUN MODEL EVIDENCE
RUN_MODEL_EV_BASE <- function(data_base, n_reps = 100){
  
  #List of model evidences
  list_log_ev = c()
  
  #REPEAT FOR REPS
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_base = MCMC_INFER_BASELINE(data_base)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_base$time_elap = time_elap
    saveRDS(mcmc_base, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base', i, '.rds' ))
    
    #MODEL EVIDENCE
    log_mod_ev = LOG_MODEL_EVIDENCE(mcmc_base$log_like_vec)
    list_log_ev = c(list_log_ev, log_mod_ev)
    print(log_mod_ev)
    
  }
  
  return(list_log_ev)
  
}

RUN_MODEL_EV_SSEB <- function(data_base, CURRENT_OUTPUT_FOLDER, n_reps = 50){
  
  #INITIALISE
  list_log_ev = c()
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_sseb = MCMC_INFER_SSEB(data_base, n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_sseb$time_elap = time_elap
    saveRDS(mcmc_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb', i, '.rds' ))
    
    #MODEL EVIDENCE
    log_mod_ev = LOG_MODEL_EVIDENCE(mcmc_sseb$log_like_vec)
    list_log_ev = c(list_log_ev, log_mod_ev)
    print(log_mod_ev)
    
  }
  
  return(list_log_ev)
}

#MODEL COMPARISON VIA POSTERIOR MODEL PROBABILITIES
#Add Do all calculations on log scale and take exp as final step 
#0.5*base_model/(0.5*base_model + (1/8)*S + (1/8)*m2 + (1/8)*m3 + (1/8)*m4)  #An 1/8 for each of 4 other models

