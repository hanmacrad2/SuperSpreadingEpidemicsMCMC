#MODEL COMPARISON 


RUN_MODEL_EVIDENCE_BASE <- function(r0, n_mcmc = 30000) {
  
  'Inference of baseline simulate data'
  
  epidemic_data = SIMULATE_EPI_BASELINE(r0)
  
  #MCMC
  mcmc_output = MCMC_INFER_BASELINE(data_baseline, n_mcmc)
  
  #Model Evidence
  log_phat = GET_LOG_MODEL_EVIDENCE_BASELINE(mcmc_output$r0_vec, epidemic_data) 
  
  return(result_row)
}