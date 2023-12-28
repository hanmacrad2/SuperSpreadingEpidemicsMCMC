#COMPUTING FUNCTIONS - INFERENCE

GET_INFER_R0_ROW <- function(r0_val, r0_vec, mcmc_output, epidemic_data){
  
  #Extract 
  num_days = length(epidemic_data)
  
  row_r0 <- data.frame(
  true_r0 = r0_val,
  mean_r0 = mean(r0_vec),
  r0_start = mcmc_output$r0_start,
  sd_r0 = sd(r0_vec),
  lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
  upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
  coverage_r0 = GET_COVERAGE(r0_val, r0_vec),
  tot_infs = sum(epidemic_data),
  end_day = epidemic_data[num_days],
  eff_size_r0 = unlist(effectiveSize(as.mcmc(r0_vec)))[[1]])
  
  row_r0$sim_data = list(epidemic_data)
    
  return(row_r0)
}

GET_INFERENCE_ALPHA <- function(true_alpha, alpha_vec){
  
  row_alpha <- data.frame(
    true_alpha = true_alpha,
    mean_alpha = mean(alpha_vec),
    sd_alpha = sd(alpha_vec),
    lower_ci_alpha = get_lower_ci(alpha_vec), # credible_intervals["lower"],
    upper_ci_alpha = get_upper_ci(alpha_vec), #credible_intervals["upper"],
    coverage_alpha = GET_COVERAGE(alpha_vec),
    esize_alpha = effectiveSize(alpha_vec)
  )
  
  return(row_alpha)
}

GET_COVERAGE <- function(param_val, mcmc_vec){
  
  #Coverage
  lower_ci = get_lower_ci(mcmc_vec)
  upper_ci = get_upper_ci(mcmc_vec)
  
  #R0
  if (param_val < lower_ci || param_val > upper_ci){
    coverage = 0
  } else {
    coverage = 1
  }
  
  return(coverage)
}

#ADD GELMAN.DIAG
#gelman.diag(as.mcmc.list(c(1,1,2,1,2,2,1,3,3,1,1,1,1,1,1,1)))
#df3 = cbind(df1, df2)


#*********************
#* . 1BASELINE MODEL
#* *******************

#MCMC FUNCTIONS
#' @export 
INFER_BASELINE  <- function(r0_val, PRIORS_USED = GET_PRIORS_USED(), n_mcmc = 40000) { #100  
  
  'Inference of baseline simulate data'
  cat(r0_val)
  epidemic_data = SIMULATE_EPI_BASELINE(r0 = r0_val)
  
  #MCMC
  mcmc_output = MCMC_INFER_BASELINE(epidemic_data, n_mcmc, PRIORS_USED = PRIORS_USED)
  r0_vec = mcmc_output$r0_vec
  
  result_row = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  
  return(result_row)
}

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
  
  result_row <- data.frame(
    true_r0 = r0_val,
    mean_r0 = mean(r0_vec),
    r0_start = mcmc_output$r0_start,
    sd_r0 = sd(r0_vec),
    lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
    upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
    coverage_r0 = GET_COVERAGE(r0_val, r0_vec),
    true_k = k_val,
    mean_k = mean(k_vec),
    sd_k = sd(k_vec),
    lower_ci_k = get_lower_ci(k_vec),
    upper_ci_k = get_upper_ci(k_vec),
    coverage_k = GET_COVERAGE(k_val, k_vec),
    tot_infs = sum(data_sse),
    end_day = data_sse[num_days],
    row.names = NULL
  )
  
  print('list data_sse')
  print(list(data_sse))
  #ADD DATA
  result_row$data_sim = list(data_sse)
  
  # Add the row to the results dataframe
  return(result_row)
  
}


#**********************
# 3. SSI MODEL
#**********************
#' @export 
INFER_SSI <- function(r0_val, k_val, PRIORS_USED, num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  data_ssi = SIMULATE_EPI_SSI(num_days = num_days, r0 = r0_val, k = k_val)
  epidemic_data = data_ssi$epidemic_data
  
  #MCMC
  mcmc_output = MCMC_INFER_SSI(epidemic_data, n_mcmc, PRIORS_USED)
  
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
    r0_start = mcmc_output$r0_start,
    sd_r0 = sd(r0_vec),
    lower_ci_r0 = get_lower_ci(r0_vec), # credible_intervals["lower"],
    upper_ci_r0 = get_upper_ci(r0_vec), #credible_intervals["upper"],
    coverage_r0 = GET_COVERAGE(r0_vec),
    true_k = k_val,
    mean_k = mean(k_vec),
    sd_k = sd(k_vec),
    lower_ci_k = get_lower_ci(k_vec),
    upper_ci_k = get_upper_ci(k_vec),
    coverage_k = GET_COVERAGE(k_vec), 
    tot_infs = sum(epidemic_data),
    end_day = epidemic_data[num_days],
    row.names = NULL
  )
  
  #Data
  result_row$data_sim <- list(epidemic_data)
  
  # Add the row to the results dataframe
  return(result_row)
  
}

#**********************************************
#* SSEB MODEL
#* ********************************************
INFER_SSEB <- function(r0_val, alpha_val, beta_val, PRIORS_USED,
                       num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0)
  data_sseb = SIMULATE_EPI_SSEB(num_days = num_days, r0 = r0_val,
                                alpha = alpha_val, beta = beta_val)
  
  # while(sum(data_sseb) < 5) {
  #   data_sseb = SIMULATE_EPI_SSEB()
  # }
  
  #MCMC
  mcmc_output = MCMC_INFER_SSEB(data_sseb, n_mcmc, PRIORS_USED)
  
  #Row result
  result_row = GET_SSEB_MCMC_ROW(r0_val, alpha_val, beta_val,
                                 mcmc_output, data_sseb, num_days)
  
  return(result_row)
}

GET_SSEB_MCMC_ROW <- function(r0_val, alpha_val, beta_val, mcmc_output,
                              data_sseb, num_days){
  
  r0_vec = mcmc_output$r0_vec
  alpha_vec = mcmc_output$alpha_vec
  gamma_vec = mcmc_output$gamma_vec
  
  result_row_r0 = GET_INFER_R0_ROW(r0_val, r0_vec, data_sseb)
  result_row_alpha = GET_INFERENCE_ALPHA(alpha_val, alpha_vec)
  
  result_row <- data.frame(
    true_alpha = alpha,
    mean_alpha = mean(alpha_vec),
    lower_ci_alpha = get_lower_ci(alpha_vec),
    upper_ci_alpha = get_upper_ci(alpha_vec),
    coverage_alpha = GET_COVERAGE(alpha_vec),
    true_gamma = gamma,
    mean_gamma = mean(gamma_vec),
    lower_ci_gamma = get_lower_ci(gamma_vec),
    upper_ci_gamma = get_upper_ci(gamma_vec),
    tot_infs = sum(data_sseb),
    end_day = data_sseb[num_days],
    row.names = NULL
  )
  
  result_row$data_sim = list(data_sseb)
  
  # Add the row to the results dataframe
  return(result_row)
  
}