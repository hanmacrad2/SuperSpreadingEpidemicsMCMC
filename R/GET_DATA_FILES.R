#GET REAL DATA FILES

GET_DATA_FILES <- function(RESULTS_FOLDER, list_models = c("baseline", "sse", "ssi", "sseb", "ssib"),
                           list_mod_ev = c("baseline", "sse", "ssi", "sseb")) {

  # List all RDS files in the directory
  files <- list.files(path = RESULTS_FOLDER, pattern = "*.rds", full.names = TRUE)
  mcmc_files <- list()
  mcmc_data_list <- list()
  
  
  #MODEL EVIDENCE
  model_evidence_files <- list()
  model_evidence_data_list <- list()
  
  # Loop through each model to filter files
  for (model in list_models) {
    mcmc_files[[model]] <- files[grepl(paste0("mcmc_", model), files)]
    model_evidence_files[[model]] <- files[grepl(paste0("mod_ev_", model), files)]
  }

  
  # Loop through each model and load the files
  for (model in list_models) {
    mcmc_data_list[[model]] <- lapply(mcmc_files[[model]], readRDS)
    model_evidence_data_list[[model]] <- lapply(model_evidence_files[[model]], readRDS)
  }
  
  
  # Loop through each model to filter files
  for (model in list_mod_ev) {
    model_evidence_files[[model]] <- files[grepl(paste0("mod_ev_", model), files)]
  }
  
  # Loop through each model and load the files
  for (model in list_mod_ev) {
    model_evidence_data_list[[model]] <- lapply(model_evidence_files[[model]], readRDS)
  }
  
  # Optionally, print the structure of the data lists to verify
  str(mcmc_data_list)
  str(model_evidence_data_list)
  
  # Access the first element for each model
  # MCMC
  mcmc_baseline <- mcmc_data_list$baseline[[1]]
  mcmc_sse <- mcmc_data_list$sse[[1]]
  mcmc_ssi <- mcmc_data_list$ssi[[1]]
  mcmc_sseb <- mcmc_data_list$sseb[[1]]
  mcmc_ssib <- mcmc_data_list$ssib[[1]]
  
  # Model Evidence
  mod_ev_baseline <- model_evidence_data_list$baseline[[1]]
  mod_ev_sse <- model_evidence_data_list$sse[[1]]
  mod_ev_ssi <- model_evidence_data_list$ssi[[1]]
  mod_ev_sseb <- model_evidence_data_list$sseb[[1]]
  #mod_ev_ssib <- model_evidence_data_list$ssib[[1]]
  
  list_mcmc = list(mcmc_baseline = mcmc_baseline, mcmc_sse = mcmc_sse, mcmc_ssi = mcmc_ssi, 
                   mcmc_sseb = mcmc_sseb, mcmc_ssib = mcmc_ssib)
  
  vec_mod_ev = c(mod_ev_baseline, mod_ev_sse, mod_ev_ssi, mod_ev_sseb) #, mod_ev_ssib)
  #list_mod_ev = list(mod_ev_baseline = mod_ev_baseline, mod_ev_sse = mod_ev_sse, mod_ev_ssi = mod_ev_ssi, 
  #                   mod_ev_sseb = mod_ev_sseb) #, mod_ev_ssib = mod_ev_ssib)
  
  return(list = list(list_mcmc = list_mcmc, vec_mod_ev = vec_mod_ev))
}