#RUN MULTIPLE MCMC
RUN_MCMC_MULTIPLE_TIMES <- function(epidemic_data, OUTER_FOLDER, run_number = 1, n_repeats = 100, n_mcmc = 50000,
                                FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                   SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #TOTAL TIME
  start_time = Sys.time()
  print(paste0('start_time:', start_time))
  
  if (FLAGS_MODELS$BASELINE){
    
    #CREATE FOLDER
    model_type = 'BASELINE'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_BASELINE(epidemic_data, n_mcmc = n_mcmc)
      
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
    
  } else if(FLAGS_MODELS$SSEB){
    
    #CREATE FOLDER
    model_type = 'SSEB'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:9){
      
      #RUN MCMC
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_SSEB(epidemic_data, n_mcmc = n_mcmc)   
      
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }

  } else if (FLAGS_MODELS$SSNB){
    
    #CREATE FOLDER
    model_type = 'SSNB'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_SSNB(epidemic_data, n_mcmc = n_mcmc)
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
    
  }
  
  #TOTAL TIME
  end_time = Sys.time()
  tot_time_elap = get_time(start_time, end_time)
  mcmc_output$tot_time_elap = tot_time_elap
  saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
}


#***************************
#* OLDER INDIVIDUAL FUNCTIONS
#***************************

#RUN SSNB MULTIPLE TIMES 
RUN_MULTIPLE_MCMC_SSNB <- function(epidemic_data, OUTER_FOLDER, model_type = 'SSNB',
                                   run_number = 1, n_reps = 100,  n_mcmc = 100000){
  
  #INITIALISE
  CURRENT_FOLDER = paste0(OUTER_FOLDER, '/', model_type, '/run_', run_number, '/')
  create_folder(CURRENT_FOLDER)
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_output = MCMC_INFER_SSNB(epidemic_data, n_mcmc = n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_output$time_elap = time_elap
    
    #SAVE
    saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', i, '.rds'))
  }
}

#RUN MULTIPLE MCMC
RUN_MULTIPLE_MCMC_SSIC <- function(epidemic_data, OUTER_FOLDER, model_type = 'SSIC',
                                   run_number = 1, n_reps = 100,  n_mcmc = 100000){
  
  #INITIALISE
  CURRENT_FOLDER = paste0(OUTER_FOLDER, '/', model_type, '/run_', run_number, '/')
  create_folder(CURRENT_FOLDER)
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_output = MCMC_INFER_SSIC(epidemic_data, n_mcmc = n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_output$time_elap = time_elap
    
    #SAVE
    saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', i, '.rds'))
  }
}

#RUN SSNB MULTIPLE TIMES
RUN_MULTIPLE_MCMC_SSIB <- function(epidemic_data, OUTER_FOLDER, model_type = 'SSIB',
                                   run_number = 1, n_reps = 100,  n_mcmc = 100000){
  
  #INITIALISE
  CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
  create_folder(CURRENT_FOLDER)
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_output = MCMC_INFER_SSIB(epidemic_data, n_mcmc = n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_output$time_elap = time_elap
    
    #SAVE
    saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i, '.rds'))
  }
}