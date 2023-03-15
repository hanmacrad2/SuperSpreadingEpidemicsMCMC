#RUN MULTIPLE MCMC
RUN_MULTIPLE_MCMC_SSEC <- function(epidemic_data, CURRENT_OUTPUT_FOLDER, model_type = 'SSEC',
                                   run_number = 1, n_reps = 100){
  
  #INITIALISE
  n_mcmc = 30000
  RESULTS_FOLDER = paste0(CURRENT_OUTPUT_FOLDER, '/', model_type, '/run_', run_number, '/')
  create_folder(RESULTS_FOLDER)
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_output = MCMC_INFER_SSEC(epidemic_data, n_mcmc = n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_output$time_elap = time_elap
    
    #SAVE
    saveRDS(mcmc_output, file = paste0(RESULTS_FOLDER, 'mcmc_', i, '.rds'))
  }
}