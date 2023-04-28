#***************************************************
#* #MODEL EVIDENCE ESTIMATE VIA HARMONIC MEAN (Kaas, Raftery, 1995)
#********************************************************
LOG_MODEL_EVIDENCE_HM <- function(loglike_vec){
  
  'Model evidence via Harmonic mean (log) (Kaas, Raftery, 1995)'
  
  #print(loglike_vec)
  loglike_vec = - loglike_vec #Harmonic mean applied to inverse likelihood
  loglike_lse = LOG_SUM_EXP(loglike_vec)
  N = length(loglike_vec)
  log_model_evidence = log(N) - loglike_lse
  
  return(log_model_evidence)
}

MODEL_EVIDENCE_HM <- function(loglike_vec){
  
  'Model evidence via Harmonic mean (Kaas, Raftery, 1995)'
  
  #print(loglike_vec)
  likelihood_vec = exp(loglike_vec);  N = length(loglike_vec)
  inner_sum = sum(1/likelihood_vec)
  harmonic_mean = 1/(inner_sum/N)
  
  return(harmonic_mean)
}

#********************************************************
#* 2. LOAD MCMC & GET MODEL EVIDENCE
#********************************************************
LOAD_MCMC_GET_MODEL_EV_HM <- function(OUTER_FOLDER, run = 2, n_repeats = 500, burn_in_pc = 0.2, BURN_IN = TRUE, #REMOVE IN FUTURE
                                      FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE,
                                                          SSNB = FALSE, SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  1. Load MCMC. 2. Get log model evidence'

  #Results
  list_log_model_ev =  rep(NA, n_repeats)

  #MODEL TYPE:
  if (FLAGS_MODELS$BASE) {
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    print(CURRENT_FOLDER)
    #model_type = 'base'
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
    #mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i))
    #mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', i ,'.rds'))
    log_like_vec = mcmc_output$log_like_vec
    
    #LOG LIKE (BURN-IN):
    if(BURN_IN){
      print(length(log_like_vec))
      log_like_vec = log_like_vec[(burn_in_pc*length(log_like_vec)):length(log_like_vec)]
      print(length(log_like_vec))
    }
    
    #LOG MODEL EVIDENCE
    list_log_model_ev[i] = LOG_MODEL_EVIDENCE_HM(log_like_vec)
  }
  
  #SAVE LOG MODEL EVIDENCE ESTIMATES
  saveRDS(list_log_model_ev, file = paste0(CURRENT_FOLDER, 'list_log_model_ev', model_type, '_', run, '.rds' ))
  
  return(list_log_model_ev) 
}

#********************************************************
#* 2.INSPECT HARMONIC MEAN ESTIMATES
#********************************************************
INSPECT_HM_MEAN_EST <- function(OUTER_FOLDER, run = 2, n_repeats = 500, burn_in_pc = 0.2, BURN_IN = TRUE, #REMOVE IN FUTURE
                                      FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE,
                                                          SSNB = FALSE, SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  1. Load MCMC. 2. Get log model evidence'
  
  #Results
  list_log_model_ev =  rep(NA, n_repeats)
  
  #MODEL TYPE:
  if (FLAGS_MODELS$BASE) {
    model_type = 'baseline'
    CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
    print(CURRENT_FOLDER)
    #model_type = 'base'
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
  #COUNTS
  count_normal = 0; count_outlier = 0
  
  #LOG MODEL EVIDENCE FOR ALL MCMC RUNS
  for (i in 1:n_repeats){
    #print(paste0('i = ', i))
    
    #MCMC OUTPUT
    mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
    log_like_vec = mcmc_output$log_like_vec
    
    #LOG LIKE (BURN-IN):
    if(BURN_IN){
      #print(length(log_like_vec))
      log_like_vec = log_like_vec[(burn_in_pc*length(log_like_vec) + 1):length(log_like_vec)]
      #print(length(log_like_vec))
    }
    
    #LOG MODEL EVIDENCE
    hm_log_model_ev = LOG_MODEL_EVIDENCE_HM(log_like_vec)
    list_log_model_ev[i] = hm_log_model_ev
    
    if(hm_log_model_ev < -92.8 & hm_log_model_ev > -92.9 & count_normal < 2){ #BASELINE:  -91.6 & hm_log_model_ev > -91.9
      print(paste0('rep = ', i))
      print('Avg case')
      print('Summary of log likelihood vector')
      print(summary(log_like_vec))
      print('')
      #print(paste0(i, ' loglike : ',  print(summary(log_like_vec)), ' log harmonic mean',  hm_log_model_ev))
      count_normal = count_normal + 1
    }
    
    if(hm_log_model_ev < -96 & count_outlier < 2){
      print(paste0('rep = ', i))
      print('Outlier, Model Evidence  < - 96 ')
      print('Summary of log likelihood vector')
      print(summary(log_like_vec))
      print('')
      #print(paste0(i, ' loglike : ',  print(summary(log_like_vec)), ' log harmonic mean',  hm_log_model_ev))
      count_outlier = count_outlier + 1
    }
  }
  
  #SAVE LOG MODEL EVIDENCE ESTIMATES
  #saveRDS(list_log_model_ev, file = paste0(CURRENT_FOLDER, 'list_log_model_ev', model_type, '_', run, '.rds' ))
  
  #print('*********'); print('**********')
  #print(list_mean_log_like)
  
  return(list_log_model_ev) 
}

