#COMPUTING FUNCTIONS - INFERENCE

#**********************
# 2. SSE MODEL
#**********************
#' @export 
INFER_SSE <- function(r0_val, k_val, num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  data_sse = SIMULATE_EPI_SSE(num_days = num_days, r0 = r0_val, k = k_val)
  
  #MCMC
  mcmc_output = MCMC_INFER_SSE(data_sse, n_mcmc)
  
  #Row result
  result_row = GET_SSE_MCMC_ROW(r0_val, k_val, mcmc_output, data_sse, num_days)
  
  return(result_row)
}

#' @export 
GET_SSE_MCMC_ROW <- function(r0_val, k_val, 
                             mcmc_output, data_sse, num_days){
  
  r0_vec = mcmc_output$sse_params_matrix[,1]
  k_vec = mcmc_output$sse_params_matrix[,2]
  
  result_row <- data.frame(
    true_r0 = r0_val,
    mean_r0 = mean(r0_vec),
    lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
    upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
    true_k = k_val,
    mean_k = mean(k_vec),
    lower_ci_k = get_lower_ci(k_vec),
    upper_ci_k = get_upper_ci(k_vec),
    tot_infs = sum(data_sse),
    end_day = data_sse[num_days],
    row.names = NULL
  )
  
  #ADD DATA
  result_row$data_sim <- list(data_sse)
  
  # Add the row to the results dataframe
  return(result_row)
  
}