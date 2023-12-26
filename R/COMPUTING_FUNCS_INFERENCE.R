#COMPUTING FUNCTIONS - INFERENCE

#**********************
# 2. SSE MODEL
#**********************
#' @export 
INFER_SSE <- function(r0_val, k_val, PRIORS_USED, num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  data_sse = SIMULATE_EPI_SSE(num_days = num_days, r0 = r0_val, k = k_val)
  
  #MCMC
  mcmc_output = MCMC_INFER_SSE(data_sse, n_mcmc, PRIORS_USED)
  
  #Row result
  result_row = GET_SSE_MCMC_ROW(r0_val, k_val, mcmc_output, data_sse, num_days)
  
  return(result_row)
}

#' @export 
GET_SSE_MCMC_ROW <- function(r0_val, k_val, 
                             mcmc_output, data_sse, num_days){
  
  r0_vec = mcmc_output$sse_params_matrix[,1]
  k_vec = mcmc_output$sse_params_matrix[,2]
  
  #Coverage
  lower_ci_r0 = get_lower_ci(r0_vec)
  upper_ci_r0 = get_upper_ci(r0_vec)
  lower_ci_k = get_lower_ci(k_vec)
  upper_ci_k = get_upper_ci(k_vec)
  
  #R0
  if (r0_val < lower_ci_r0 || r0_val > upper_ci_r0){
    coverage_in_r0 = 0
  } else {
    coverage_in_r0 = 1
  }
  
  #R0
  if (k_val < lower_ci_k || k_val > upper_ci_k){
    coverage_in_k = 0
  } else {
    coverage_in_k = 1
  }
  
  result_row <- data.frame(
    true_r0 = r0_val,
    mean_r0 = mean(r0_vec),
    r0_start = mcmc_output$r0_start,
    sd_r0 = sd(r0_vec),
    lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
    upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
    coverage_r0 = coverage_in_r0,
    true_k = k_val,
    mean_k = mean(k_vec),
    sd_k = sd(k_vec),
    lower_ci_k = get_lower_ci(k_vec),
    upper_ci_k = get_upper_ci(k_vec),
    coverage_k = coverage_in_k,
    tot_infs = sum(data_sse),
    end_day = data_sse[num_days],
    data_sim = list(data_sse),
    row.names = NULL
  )
  
  #ADD DATA
  #result_row$data_sim <- list(data_sse)
  
  # Add the row to the results dataframe
  return(result_row)
  
}


#**********************
# 3. SSI MODEL
#**********************
#' @export 
INFER_SSI <- function(r0_val, k_val, num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  data_ssi = SIMULATE_EPI_SSI(num_days = num_days, r0 = r0_val, k = k_val)
  epidemic_data = data_ssi$epidemic_data
  
  #MCMC
  mcmc_output = MCMC_INFER_SSI(epidemic_data, n_mcmc)
  
  #Row result
  result_row = GET_SSI_MCMC_ROW(r0_val, k_val, mcmc_output, epidemic_data, num_days)
  
  return(result_row)
}

#' @export 
GET_SSI_MCMC_ROW <- function(r0_val, k_val, 
                             mcmc_output, epidemic_data, num_days){
  
  r0_vec = mcmc_output$ssi_params_matrix[,1]
  k_vec = mcmc_output$ssi_params_matrix[,2]
  
  result_row <- data.frame(
    true_r0 = r0_val,
    mean_r0 = mean(r0_vec),
    lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
    upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
    true_k = k_val,
    mean_k = mean(k_vec),
    lower_ci_k = get_lower_ci(k_vec),
    upper_ci_k = get_upper_ci(k_vec),
    tot_infs = sum(epidemic_data),
    end_day = epidemic_data[num_days],
    row.names = NULL
  )
  
  #ADD DATA
  result_row$data_sim <- list(epidemic_data)
  
  # Add the row to the results dataframe
  return(result_row)
  
}