#******************
#MULTIPLE RJMCMC
#****************
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