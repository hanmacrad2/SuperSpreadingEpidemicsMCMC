#RUN MULTIPLE MCMC

#RUN SSNB MULTIPLE TIMES 
RUN_MULTIPLE_MCMC_SSNB <- function(epidemic_data, OUTER_FOLDER, model_type = 'SSNB',
                                   run_number = 1, n_reps = 100,  n_mcmc = 100000){
  
  #INITIALISE
  RESULTS_FOLDER = paste0(OUTER_FOLDER, '/', model_type, '/run_', run_number, '/')
  create_folder(RESULTS_FOLDER)
  
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
    saveRDS(mcmc_output, file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
  }
}

#RUN MULTIPLE MCMC
RUN_MULTIPLE_MCMC_SSIC <- function(epidemic_data, OUTER_FOLDER, model_type = 'SSIC',
                                   run_number = 1, n_reps = 100,  n_mcmc = 100000){
  
  #INITIALISE
  RESULTS_FOLDER = paste0(OUTER_FOLDER, '/', model_type, '/run_', run_number, '/')
  create_folder(RESULTS_FOLDER)
  
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
    saveRDS(mcmc_output, file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
  }
}

#RUN SSNB MULTIPLE TIMES
RUN_MULTIPLE_MCMC_SSIB <- function(epidemic_data, OUTER_FOLDER, model_type = 'SSIB',
                                   run_number = 1, n_reps = 100,  n_mcmc = 100000){
  
  #INITIALISE
  RESULTS_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
  create_folder(RESULTS_FOLDER)
  
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
    saveRDS(mcmc_output, file = paste0(RESULTS_FOLDER, 'mcmc_', tolower(model_type), '_', i, '.rds'))
  }
}