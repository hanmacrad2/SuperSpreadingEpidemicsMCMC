#RUN SSE PARALLEL
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
library(MASS)
library(extraDistr)
library(RChronoModel)
library(dplyr)
library(parallel)

#SCRIPT
cat('RUNNING SCRIPT \n')

# Define a custom function to print the iteration number
print_iteration <- function(iteration) {
  cat("Processing iteration", iteration, "\n")
}

#**************
#DATA
#**************

#SARS 03
data_type = 'sars_canada_03'
file_data = 'epi_data_sars_canada_03.rds'
epidemic_data = readRDS(file_data)

#FOLDER SAVE
#FOLDER = '~/GitHub/computing/REAL_DATA/HONG_KONG/MCMC/WAVE_2/'
#FOLDER = './MCMC/' #Same thing
#create_folder(FOLDER)

RUN_MODELS_REAL_DATA <- function(iteration, epidemic_data, n_mcmc = 62500){
  
  # Initialize the vector for model evidence
  vec_mod_ev = c()
  
  #***************************
  #1. BASELINE MODEL
  if(iteration == 1){
    
    mcmc_baseline = MCMC_INFER_BASELINE(epidemic_data, n_mcmc)
    time_stp = GET_CURRENT_TIME_STAMP()
    file_name = paste0('mcmc_baseline_', time_stp, '_', data_type, '.rds')
    saveRDS(mcmc_baseline, file = file_name)
    
    #MODEL EVIDENCE 
    file_name = paste0('mod_ev_baseline_', time_stp, '_', data_type, '.rds')
    mod_ev_base = GET_LOG_MODEL_EVIDENCE_BASELINE(mcmc_baseline$r0_vec, epidemic_data)
    mod_ev_base
    saveRDS(mod_ev_base, file = file_name)
    vec_mod_ev = c(vec_mod_ev, mod_ev_base)
    
  } else if (iteration == 2){
    
    #*****************************
    #2. SSE MODEL
    mcmc_sse = MCMC_INFER_SSE(epidemic_data, n_mcmc)
    time_stp = GET_CURRENT_TIME_STAMP()
    file_name = paste0('mcmc_sse_', time_stp, '_', data_type, '.rds')
    saveRDS(mcmc_sse, file = file_name)
    
    #MODEL EVIDENCE
    mcmc_samples =  mcmc_sse$sse_params_matrix 
    FLAGS_MODELS = GET_FLAGS_MODELS(SSE = TRUE) 
    mod_ev_sse = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data, FLAGS_MODELS = FLAGS_MODELS) 
    mod_ev_sse
    file_name = paste0('mod_ev_sse_', time_stp, '_', data_type, '.rds')
    saveRDS(mod_ev_sse, file = file_name)
    vec_mod_ev = c(vec_mod_ev, mod_ev_sse)
    
    
  } else if (iteration == 3){
    
    #*****************************
    #3. SSI MODEL
    mcmc_ssi = MCMC_INFER_SSI(epidemic_data, n_mcmc)
    
    time_stp = GET_CURRENT_TIME_STAMP()
    file_name = paste0('mcmc_ssi_', time_stp, '_', data_type, '.rds')
    saveRDS(mcmc_ssi, file = file_name)
    
    #MODEL EVIDENCE
    mod_ev_ssi = GET_LOG_MODEL_EVIDENCE_SSI(mcmc_ssi, epidemic_data) 
    mod_ev_ssi
    file_name = paste0('mod_ev_ssi_', time_stp, '_', data_type, '.rds')
    saveRDS(mod_ev_ssi, file = file_name)
    vec_mod_ev = c(vec_mod_ev, mod_ev_ssi)
    
  } else if (iteration == 4){
    
    #*****************************
    #4. SSEB
    mcmc_sseb = MCMC_INFER_SSEB(epidemic_data, n_mcmc)
    
    time_stp = GET_CURRENT_TIME_STAMP()
    file_name = paste0('mcmc_sseb_', time_stp, '_', data_type, '.rds')
    saveRDS(mcmc_sseb, file = file_name)
    
    #MODEL EVIDENCE
    mcmc_samples =  matrix(c(mcmc_sseb$alpha_vec, mcmc_sseb$r0_vec,
                             mcmc_sseb$beta_vec), ncol = 3)
    
    FLAGS_MODELS = GET_FLAGS_MODELS(SSEB = TRUE) 
    
    mod_ev_sseb = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data,
                                         FLAGS_MODELS = FLAGS_MODELS) 
    mod_ev_sseb
    file_name = paste0('mod_ev_sseb_', time_stp, '_', data_type, '.rds')
    saveRDS(mod_ev_sseb, file = file_name)
    vec_mod_ev = c(vec_mod_ev, mod_ev_sseb)
    
  } else if (iteration == 5){
    
    #***********
    #* 5. SSIB
    mcmc_ssib = MCMC_INFER_SSIB(epidemic_data, n_mcmc)
    
    time_stp = GET_CURRENT_TIME_STAMP()
    file_name = paste0('mcmc_ssib_', time_stp, '_', data_type, '.rds')
    saveRDS(mcmc_ssib, file = file_name)
    
    #MODEL EVIDENCE
    mod_ev_ssib = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_ssib, epidemic_data)
    
    mod_ev_ssib
    file_name = paste0('mod_ev_ssib_', time_stp, '_', data_type, '.rds')
    saveRDS(mod_ev_sseb, file = file_name)
    vec_mod_ev = c(vec_mod_ev, mod_ev_ssib)
    
  }
  
  return(vec_mod_ev) 
}

#APPLY
num_runs = 5
# Use mclapply to run the INFER_SSIB function in parallel
vec_mod_ev_list  <- mclapply(1:num_runs, function(iteration) {
  tryCatch({
    print_iteration(iteration)  # Print the iteration number
    RUN_MODELS_REAL_DATA(iteration, epidemic_data)
  }, error = function(e) {
    cat("Error in GET_MODEL_EVIDENCE for run =", iteration, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    NULL
  })
})

vec_mod_ev = c()
# Concatenate the model evidences into a single vector
for (evidence in vec_mod_ev_list) {
  if (!is.null(evidence)) {
    vec_mod_ev <- c(vec_mod_ev, unlist(evidence))
  }
}


#*********************************
#*
# RUN POSTERIOR MODEL PROBABILITES
#*
#**********************************

time_stp = GET_CURRENT_TIME_STAMP()
file_name = paste0('vec_mod_ev', time_stp, '_', data_type, '.rds')
saveRDS(vec_mod_ev, file = file_name)

post_probs = GET_MODELS_POSTERIOR_PROBS(vec_mod_ev)
file_name = paste0('post_probs_', time_stp, '_', data_type, '.rds')
saveRDS(post_probs, file = file_name)
