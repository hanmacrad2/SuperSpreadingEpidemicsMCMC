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
RUN_MODEL_EV_BASE <- function(epidemic_data, CURRENT_OUTPUT_FOLDER, n_reps = 30){
  
  #List of model evidences
  list_log_ev = c()
  
  #REPEAT FOR REPS
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_base = MCMC_INFER_BASELINE(epidemic_data)
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

#MODEL EVIDENCE - SSEB
RUN_MODEL_EV_SSEB <- function(epidemic_data, CURRENT_OUTPUT_FOLDER, n_reps = 30){
  
  #INITIALISE
  n_mcmc = 30000
  list_log_ev = c()
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_sseb = MCMC_INFER_SSEB(epidemic_data, n_mcmc)
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

#MULTIPLE RJMCMC
RUN_RJMCMC_MULT <- function(epidemic_data, CURRENT_OUTPUT_FOLDER, n_reps = 10){
  
  #INITIALISE
  n_mcmc = 30000
  list_bfs = c(); list_bc0 = c()
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    rj_output = RJMCMC_BASE_SSEB(epidemic_data, n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    rj_output$time_elap = time_elap
    saveRDS(rj_output, file = paste0(CURRENT_OUTPUT_FOLDER, '/rjmcmc_', i, '.rds' ))
    
    #MODEL EVIDENCE
    list_bfs = c(list_bfs, rj_output$bayes_factor)
    print(list_bfs)
    
    list_bc0 = c(list_bc0, rj_output$beta_pc0)
    print(list_bc0)
    
  }
  
  rjmcmc_out = list(list_bfs = list_bfs, list_bc0 = list_bc0)
  
  return(rjmcmc_out)
}

#MODEL COMPARISON VIA POSTERIOR MODEL PROBABILITIES
#Add Do all calculations on log scale and take exp as final step 
#(Not one 8th)
#0.5*base_model/(0.5*base_model + (1/8)*S + (1/8)*m2 + (1/8)*m3 + (1/8)*m4)  #An 1/8 for each of 4 other models

#Check
list1 = c()
for(i in 1:10){
  
  list1 = c(list1, i**2)
  #print(list1)

}
list1

list2 = list(list1 = list1)
