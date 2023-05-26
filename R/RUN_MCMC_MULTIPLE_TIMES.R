#RUN MULTIPLE MCMC
RUN_MCMC_MULTIPLE_TIMES <- function(epidemic_data, OUTER_FOLDER, 
                                    run_number = 1, n_repeats = 100, n_mcmc = 50000,
                                    PRIORS_USED = list(EXP_K = TRUE, GAMMA_K = FALSE),
                                FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                  SSIR = FALSE, SSIB = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #TOTAL TIME
  start_time = Sys.time()
  print(paste0('start_time:', start_time))
  print(FLAGS_MODELS)
  
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
      print(mean(mcmc_output$r0_vec))
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
    
  } else if(FLAGS_MODELS$SSEB){
    
    #CREATE FOLDER
    model_type = 'SSEB'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_SSEB(epidemic_data, n_mcmc = n_mcmc)   
      print(mean(mcmc_output$alpha_vec))
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }

  } else if (FLAGS_MODELS$SSNB){
    
    #CREATE FOLDER
    model_type = 'SSNB'; print(model_type); print(PRIORS_USED)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC 
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_SSNB(epidemic_data, n_mcmc = n_mcmc, PRIORS_USED = PRIORS_USED)
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
    
  } else if (FLAGS_MODELS$SSIR) {
    
    #CREATE FOLDER
    model_type = 'SSIR'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_SSIR(epidemic_data, n_mcmc = n_mcmc)
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
  } else if (FLAGS_MODELS$SSIB) {
    
    #CREATE FOLDER
    model_type = 'SSIB'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      mcmc_output = MCMC_INFER_SSIB(epidemic_data, n_mcmc = n_mcmc)
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

#***********
#* RUN_MCMC_MULTIPLE_DATASETS
#* **************************
RUN_MCMC_MULTIPLE_DATASETS <- function(matrix_datasets, OUTER_FOLDER, run_number = 1, n_mcmc = 20000,
                                       PRIORS_USED = list(EXP_K = TRUE, GAMMA_K = FALSE),
                                    FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                        SSIB = FALSE, SSIR = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #TOTAL TIME
  n_repeats = dim(matrix_datasets)[1]
  start_time = Sys.time()
  print(paste0('start_time:', start_time))
  
  if (FLAGS_MODELS$BASELINE){
    
    #CREATE FOLDER
    model_type = 'BASELINE'; print(model_type)
    
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type , '/') 
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      epidemic_data = matrix_datasets[i,]
      mcmc_output = MCMC_INFER_BASELINE(epidemic_data, n_mcmc = n_mcmc)
      print(mean(mcmc_output$r0_vec))
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
    
  }
  
  if(FLAGS_MODELS$SSEB){
    
    #CREATE FOLDER
    model_type = 'SSEBX'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      epidemic_data = matrix_datasets[i,]
      mcmc_output = MCMC_INFER_SSEB(epidemic_data, n_mcmc = n_mcmc)   
      #print(mean(mcmc_output$alpha_vec))
      #SAVE
      saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))
    }
    
  }
  
  if (FLAGS_MODELS$SSNB){
    
    #CREATE FOLDER
    model_type = 'SSNB'; print(model_type)
    CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/')
    create_folder(CURRENT_FOLDER)
    print(paste0('CURRENT_FOLDER = ', CURRENT_FOLDER))
    
    for (i in 1:n_repeats){
      
      #RUN MCMC
      print(paste0('i = ', i))
      epidemic_data = matrix_datasets[i,]
      mcmc_output = MCMC_INFER_SSNB(epidemic_data, n_mcmc = n_mcmc, PRIORS_USED = PRIORS_USED)
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


