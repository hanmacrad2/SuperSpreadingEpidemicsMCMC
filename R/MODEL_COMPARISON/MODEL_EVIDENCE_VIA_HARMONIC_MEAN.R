#***************************************************
#* #MODEL EVIDENCE ESTIMATE VIA HARMONIC MEAN (Kaas, Raftery, 1995)
#********************************************************
LOG_MODEL_EVIDENCE_HM <- function(loglike_vec){
  
  'Model evidence via Harmonic mean (log) (Kaas, Raftery, 1995)'
  
  loglike_vec = - loglike_vec #Harmonic mean applied to inverse likelihood
  loglike_lse = LOG_SUM_EXP(loglike_vec)
  N = length(loglike_vec)
  log_model_evidence = log(N) - loglike_lse
  
  return(log_model_evidence)
}

MODEL_EVIDENCE_HM <- function(loglike_vec){
  
  'Model evidence via Harmonic mean (Kaas, Raftery, 1995)'
  
  likelihood_vec = exp(loglike_vec);  N = length(loglike_vec)
  inner_sum = sum(1/likelihood_vec)
  harmonic_mean = 1/(inner_sum/N)
  
  return(harmonic_mean)
}

#********************************************************
#* 2. LOAD MCMC & GET MODEL EVIDENCE
#********************************************************
LOAD_MCMC_GET_MODEL_EV_HM <- function(OUTER_FOLDER, run = 1, n_repeats = 100, burn_in_pc = 0.2, BURN_IN = TRUE, #REMOVE IN FUTURE
                                      FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                          SSNB = FALSE, SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  1. Load MCMC. 2. Get log model evidence'

  #Results
  list_log_model_ev =  rep(NA, n_repeats)
  list_mean_log_like =  rep(NA, n_repeats)

  #MODEL TYPE:
  if (FLAGS_MODELS$BASE) {
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  } else if (FLAGS_MODELS$SSEB)  {
    model_type = 'sseb'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  } else if  (FLAGS_MODELS$SSIB)  {
    model_type = 'ssib'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  } else if (FLAGS_MODELS$SSNB)  {
    model_type = 'ssnb'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    print(CURRENT_FOLDER)
  }
  
  #LOG MODEL EVIDENCE FOR ALL MCMC RUNS
  for (i in 1:n_repeats){
    print(paste0('i = ', i))
    
    #MCMC OUTPUT
    mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
    log_like_vec = mcmc_output$log_like_vec
    
    #LOG LIKE (BURN-IN):
    if(BURN_IN){
      print(length(log_like_vec))
      log_like_vec = log_like_vec[(burn_in_pc*length(log_like_vec)):length(log_like_vec)]
      print(length(log_like_vec))
    }
    
    #LOG MODEL EVIDENCE
    list_log_model_ev[i] = LOG_MODEL_EVIDENCE_HM(log_like_vec)
    #list_mean_log_like[i] = mean(log_like_vec)
    print(paste0('loglike : ',  print(summary(log_like_vec)), ' log harmonic mean',   list_log_model_ev[i]))
    #print(list_log_model_ev)
  }
  
  #SAVE LOG MODEL EVIDENCE ESTIMATES
  saveRDS(list_log_model_ev, file = paste0(CURRENT_FOLDER, 'list_log_model_ev', model_type, '_', run, '.rds' ))
  
  #print('*********'); print('**********')
  #print(list_mean_log_like)
  
  return(list_log_model_ev) 
}

#********************************************************
#* OLDER FUNCTIONS
#********************************************************
RUN_MODEL_EV_HM_BASE <- function(epidemic_data, CURRENT_OUTPUT_FOLDER, n_reps = 30){
  
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
    log_mod_ev = LOG_MODEL_EVIDENCE_HM(mcmc_base$log_like_vec)
    list_log_ev = c(list_log_ev, log_mod_ev)
    print(log_mod_ev)
    
  }
  
  return(list_log_ev)
  
}

#MODEL EVIDENCE - SSEB
RUN_MODEL_EV_HM_SSEB <- function(epidemic_data, CURRENT_OUTPUT_FOLDER, n_reps = 30){
  
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
    log_mod_ev = LOG_MODEL_EVIDENCE_HM(mcmc_sseb$log_like_vec)
    list_log_ev = c(list_log_ev, log_mod_ev)
    print(log_mod_ev)
    
  }
  
  return(list_log_ev)
}