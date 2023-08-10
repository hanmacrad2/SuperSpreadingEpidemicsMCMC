#!/software/easybuild/software/R/4.2.1-foss-2022a/bin/R
library(SuperSpreadingEpidemicsMCMC)
library(RChronoModel)
library(dplyr)

#PARAMS
n_mcmc = 30000

#SCRIPT
print('RUNNING SCRIPT')
args <- commandArgs(trailingOnly = TRUE)
array_index <- as.integer(args[1])
print(paste0('array_index: ', array_index))
set.seed(array_index)

#*********************
#FUNCTIONS
#*********************
GET_CURRENT_TIME_STAMP <- function(){
  
  current_time <- Sys.time()
  time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  
  return(time_string)
}

GET_FOLDER_TIME_STAMP <- function(folder_type, array_index){
  
  current_time <- Sys.time()
  time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  CURRENT_FOLDER <- paste0(folder_type, '_', array_index, '_', time_string, '/')
  print(CURRENT_FOLDER)
  
  return(CURRENT_FOLDER)
}

#MCMC
GET_MCMC_ROW <- function(true_r0, mcmc_output){
  
  # Create a row for the dataframe
  r0_vec = mcmc_output$r0_vec
  
  result_row <- data.frame(
    true_r0 = true_r0,
    lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
    upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
    row.names = NULL
  )
  
  #Make numeric
  result_row$true_r0 <- as.numeric(result_row$true_r0)
  result_row$lower_ci_r0 <- as.numeric(result_row$lower_ci_r0)
  result_row$upper_ci_r0 <- as.numeric(result_row$upper_ci_r0)
  
  # Add the row to the results dataframe
  return(result_row)
  
}

#- Create a data frame 
RUN_INFERENCE_BASE <- function(CURRENT_FOLDER, num_runs = 101, n_mcmc = 30000) {
  
  df_results = data.frame()
  list_r0 = seq(from = 0.8, to = 2.0, length.out = num_runs)
  
  for (i in 1:length(list_r0)){
    
    r0 = list_r0[i]
    data_baseline = SIMULATE_EPI_BASELINE(r0)
    data_file <- paste("data_baseline_", r0, "_", i, ".rds", sep = "")
    
    #MCMC
    mcmc_output = MCMC_INFER_BASELINE(data_baseline, n_mcmc = n_mcmc)
    
    #SAVE
    mcmc_file <- paste("mcmc_baseline_", r0, "_", i, ".rds", sep = "")
    
    #Row result
    result_row = GET_MCMC_ROW(r0, mcmc_output)
    df_results <- bind_rows(df_results, result_row)
    
    #SAVE
    saveRDS(data_baseline, file = paste0(CURRENT_FOLDER, data_file))
    saveRDS(mcmc_output, file = paste0(CURRENT_FOLDER, mcmc_file))
    
  }
  
  #SAVE DATAFRAME 
  time_stp = GET_CURRENT_TIME_STAMP()
  saveRDS(df_results, file = paste0('base_mcmc_', time_stp, '.rds'))
}


#************************
#* IMPLEMENT
#***********************

#FOLDER
folder_type = 'data_base'
CURRENT_FOLDER = GET_FOLDER_TIME_STAMP(folder_type, array_index)
create_folder(CURRENT_FOLDER)

#RUN_INFERENCE_BASE
RUN_INFERENCE_BASE(CURRENT_FOLDER)


