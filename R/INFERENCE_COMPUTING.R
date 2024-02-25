#******************************************
#COMPUTING FUNCTIONS - INFERENCE
#******************************************

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
  row_r0$bias_r0 = row_r0$mean_r0 - row_r0$true_r0
  row_r0$MAE_r0 = abs(row_r0$mean_r0 - row_r0$true_r0)
    
  return(row_r0)
}

GET_INFER_K_ROW <- function(k_val, k_vec){
  
  row_k <- data.frame(
    true_k = k_val,
    mean_k = mean(k_vec),
    sd_k = sd(k_vec),
    lower_ci_k = get_lower_ci(k_vec),
    upper_ci_k = get_upper_ci(k_vec),
    coverage_k = GET_COVERAGE(k_val, k_vec),
    eff_size_k = unlist(effectiveSize(as.mcmc(k_vec)))[[1]])
  
  #Bias
  row_k$bias_k = row_k$mean_k - row_k$true_k
  row_k$MAE_k = abs(row_k$mean_k - row_k$true_k)
  
  return(row_k)
}

#PARAMTER INFERENCE
GET_PARAM_INFERENCE <- function(true_val, param_vec, FLAG_PARAM){
  
  # Get names of true parameters
  param_name <- names(FLAG_PARAM)[unlist(FLAG_PARAM)]
  
  # Create a function to generate column names
  generate_column_names <- function(prefix, param_name) {
    
    return(paste(prefix, param_name, sep = "_"))
  }
  
  #Names in dataframe - for different statistics
  col_names <- c(
    generate_column_names("true", param_name),
    generate_column_names("mean", param_name),
    generate_column_names("sd", param_name),
    generate_column_names("lower_ci", param_name),
    generate_column_names("upper_ci", param_name),
    generate_column_names("coverage", param_name),
    generate_column_names("eff_size", param_name)
  )
  
  #Values; dataframe
  values_list <- list(
    true_val,
    mean(param_vec),
    sd(param_vec),
    get_lower_ci(param_vec),
    get_upper_ci(param_vec),
    GET_COVERAGE(true_val, param_vec),
    effectiveSize(param_vec)
  )
  
  row_result <- data.frame(matrix(unlist(values_list), nrow = 1))
  names(row_result) <- col_names
  
  row_result[paste0('bias_', param_name)] = row_result[,paste0('mean_', param_name)] - row_result[paste0('true_', param_name)]
  row_result[paste0('MAE_', param_name)] = abs(row_result[,paste0('mean_', param_name)] - row_result[paste0('true_', param_name)])
  
  return(row_result)
}

GET_FLAG_PARAM <- function(){
  
  FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE,
                    alpha = FALSE, a = FALSE, 
                    beta = FALSE, b = FALSE, c = FALSE)
  
  return(FLAG_PARAM)
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

#***************************************************************
#* 1. BASELINE MODEL
#**************************************************************

#' @export 
INFER_BASELINE  <- function(r0_val, n_mcmc, PRIORS_USED,
                            STORE_MCMC = TRUE, lt_val = 20) { #100  
  
  'Inference of baseline simulate data'
  cat(r0_val)
  
  #DATA
  epidemic_data = SIMULATE_EPI_BASELINE(r0 = r0_val)
  
  while(sum(epidemic_data) < lt_val){
    epidemic_data = SIMULATE_EPI_BASELINE(r0 = r0_val)
  }
  
  #MCMC
  mcmc_output = MCMC_INFER_BASELINE(epidemic_data, n_mcmc, PRIORS_USED = PRIORS_USED)
  r0_vec = mcmc_output$r0_vec
  
  result_row = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  result_row$accept_rate = mcmc_output$accept_rate
    
  if(STORE_MCMC){
    result_row$r0_mcmc = list(r0_vec)
  }
  
  return(result_row)
}

#*******************************************
# 2. SSE MODEL
#*******************************************
#' @export 
INFER_SSE <- function(r0_val, k_val, n_mcmc, lt_val = 20, STORE_MCMC = TRUE) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  
  epidemic_data = SIMULATE_EPI_SSE(r0 = r0_val, k = k_val)
  
  while(sum(epidemic_data) < lt_val){
    epidemic_data = SIMULATE_EPI_SSE(r0 = r0_val, k = k_val)
  }
  
  #MCMC
  mcmc_output = MCMC_INFER_SSE(epidemic_data, n_mcmc)
  r0_vec = mcmc_output$sse_params_matrix[,1]
  k_vec = mcmc_output$sse_params_matrix[,2]
  
  #Row result
  r0_row = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  k_row = GET_INFER_K_ROW(k_val, k_vec)
  
  if (STORE_MCMC){
    r0_row$r0_mcmc = list(r0_vec)
    k_row$k_mcmc = list(k_vec) 
  }
  
  result_row <- cbind(r0_row, k_row)
  result_row$accept_rate = mcmc_output$accept_rate
  
  return(result_row)
}

#*******************************************
# 3. SSI MODEL
#*******************************************
#' @export 
INFER_SSI <- function(r0_val, k_val, n_mcmc, STORE_MCMC = TRUE){ #{40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  data_ssi = SIMULATE_EPI_SSI(r0 = r0_val, k = k_val)
  epidemic_data = data_ssi$epidemic_data
  
  while(sum(epidemic_data)< 30){
    data_ssi = SIMULATE_EPI_SSI(r0 = r0_val, k = k_val)
    epidemic_data = data_ssi$epidemic_data
  }
  
  #MCMC
  mcmc_output = MCMC_INFER_SSI(epidemic_data, n_mcmc)
  r0_vec = mcmc_output$ssi_params_matrix[,1]
  k_vec = mcmc_output$ssi_params_matrix[,2]

  #Row result
  r0_row = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  k_row = GET_INFER_K_ROW(k_val, k_vec)
  
  if (STORE_MCMC){
    r0_row$r0_mcmc = list(r0_vec)
    k_row$k_mcmc = list(k_vec) 
  }

  result_row <- cbind(r0_row, k_row)
  result_row$accept_rate = mcmc_output$accept_rate
  
  return(result_row)
}


#**********************************************
#* SSEB MODEL
#* ********************************************
INFER_SSEB <- function(r0_val, alpha_val, beta_val, n_mcmc, STORE_MCMC = TRUE) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  epidemic_data = SIMULATE_EPI_SSEB(r0 = r0_val,
                                alpha = alpha_val, beta = beta_val)
  
  while(sum(epidemic_data) < 30){
    epidemic_data = SIMULATE_EPI_SSEB(r0 = r0_val,
                                      alpha = alpha_val, beta = beta_val)
  }
  
  #MCMC
  mcmc_output = MCMC_INFER_SSEB(epidemic_data, n_mcmc)
  
  #Row result - parameters
  row_r0 = GET_INFER_R0_ROW(r0_val, mcmc_output$r0_vec, mcmc_output, epidemic_data)
  
  #alpha
  FLAG_PARAM = GET_FLAG_PARAM()
  FLAG_PARAM$alpha = TRUE
  row_alpha = GET_PARAM_INFERENCE(alpha_val, mcmc_output$alpha_vec, FLAG_PARAM)
  
  #beta
  FLAG_PARAM = GET_FLAG_PARAM()
  FLAG_PARAM$beta = TRUE
  row_beta = GET_PARAM_INFERENCE(beta_val, mcmc_output$beta_vec, FLAG_PARAM)
  
  #COMBINE
  result_row = cbind(row_r0, row_alpha, row_beta)
  
  #ACCEPTANCE RATES
  result_row$accept_rate_r0 = round(mcmc_output$list_accept_rates$accept_rate_r0, 3)
  result_row$accept_rate_alpha = round(mcmc_output$list_accept_rates$accept_rate_alpha, 3)
  result_row$accept_rate_beta = round(mcmc_output$list_accept_rates$accept_rate_beta, 3)
  
  if (STORE_MCMC){
    result_row$r0_mcmc = list(mcmc_output$r0_vec)
    result_row$alpha_mcmc = list(mcmc_output$alpha_vec)
    result_row$beta_mcmc = list(mcmc_output$beta_vec)
  }
  
  
  return(result_row)
}


#**********************************************
#* SSIB 
#************************************************

INFER_SSIB <- function(r0_val, a_val, b_val, n_mcmc, STORE_MCMC = TRUE) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  epidemic_data = SIMULATE_EPI_SSIB(r0 = r0_val,
                                    a = a_val, b = b_val) 
  
  while(sum(epidemic_data) < 30){
    epidemic_data = SIMULATE_EPI_SSIB(r0 = r0_val, a = a_val, b = b_val)
  }
  
  #SSIB; pass in data
  # list_ssib_data = SIMULATE_EPI_SSIB_LIST(num_days = num_days, r0 = r0_val,
  #                                    a = a_val, b = b_val)
  # 
  # data = list(non_ss = list_ssib_data$non_ss, ss = list_ssib_data$ss)
  # epidemic_data = as.vector(unlist(list_ssib_data$total_infections))
  
  #MCMC  
  #mcmc_output = MCMC_INFER_SSIB_JOINT_II(epidemic_data, list_ssib_data, n_mcmc)
  mcmc_output = MCMC_INFER_SSIB(epidemic_data, n_mcmc)
  
  #PARAMS
  ssib_params_matrix = mcmc_output$ssib_params_matrix
  r0_vec = ssib_params_matrix[,1]
  a_vec = ssib_params_matrix[,2]
  b_vec = ssib_params_matrix[,3]
  
  #Row result - parameters
  row_r0 = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  #row_r0$mcmc_r0 = list(mcmc_output$r0_vec)
  
  #alpha
  FLAG_PARAM = GET_FLAG_PARAM()
  FLAG_PARAM$a = TRUE
  row_a = GET_PARAM_INFERENCE(a_val, a_vec, FLAG_PARAM)
  
  #beta
  FLAG_PARAM = GET_FLAG_PARAM()
  FLAG_PARAM$b = TRUE
  row_b = GET_PARAM_INFERENCE(b_val, b_vec, FLAG_PARAM)
  
  #browser()
  result_row = cbind(row_r0, row_a, row_b)
  
  #SS Data
  result_row$ss_data = list(mcmc_output$ss_inf)
  result_row$ns_data = list(mcmc_output$ns_inf)
  result_row$accept_rate = mcmc_output$accept_rate
  result_row$accept_da = mcmc_output$accept_da
  #result_row$accept_da = list(mcmc_output$vec_accept_da)
  
  if (STORE_MCMC){
    result_row$r0_mcmc = list(r0_vec)
    result_row$a_mcmc = list(a_vec)
    result_row$b_mcmc = list(b_vec)
  }
  
  return(result_row)
}

#*************************************************
#*
#* PLOT PERFORMANCE & INFERENCE RESULTS
#* 
#* ************************************************
SIM_PERFORMANCE_R0 <- function(df_results){
  
  #Bias, MAE
  df_results$MAE = abs(df_results$mean_r0 - df_results$true_r0)
  df_results$bias = df_results$mean_r0 - df_results$true_r0
  num_runs = length(df_results$true_r0)
  
  #Results
  print(paste0('mean bias: ', round(mean(df_results$bias), 3)))
  print(paste0('MAE: ', round(mean(df_results$MAE), 3)))
  print(paste0('mean sd: ', round(mean(df_results$sd_r0), 3)))
  print(paste0('coverage: ', sum(df_results$coverage_r0)))
  print(paste0('% coverage: ', sum(df_results$coverage_r0)/num_runs))
  
  #return(df_results)
  
}

#*********************************
#*
#* PERFORMANCE METRICS
#* 
#* **********************************
SIM_PERFORMANCE <- function(df_results, FLAG_PARAM = GET_PARAM(r0 = TRUE), SSEB = FALSE){
    
    #Param
    param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
    print(paste0(param, ' Performance metrics'))
    num_runs = length(df_results$true_r0)
    
    #Bias, MAE, coverage 
    mean_est = mean(unlist(df_results[paste0('mean_', param)]))
    mean_eff = mean(unlist(df_results[paste0('eff_size_', param)]))
    lower_ci_mean = mean(unlist(df_results[paste0('lower_ci_', param)]))
    upper_ci_mean = mean(unlist(df_results[paste0('upper_ci_', param)]))
    MAE = unlist(as.vector(abs(df_results[paste0('true_', param)] - df_results[paste0('mean_', param)])))
    bias = unlist(as.vector(df_results[paste0('true_', param)] - df_results[paste0('mean_', param)]))
    sd = unlist(as.vector(df_results[paste0('sd_', param)]))
    coverage = unlist(as.vector(df_results[paste0('coverage_', param)]))
    coverage_pc = sum(coverage)/num_runs
    
    if(SSEB){
      mean_accept_rate = mean(unlist(df_results[paste0('accept_rate_', param)]))
      print(paste0('accept rate ', param, ':', round(mean_accept_rate, 3))) 
    } else {
      accept_rate = mean(unlist(as.vector(df_results['accept_rate'])))
      print(paste0('accept rate: ', round(accept_rate, 3))) 
    }
   
    #Results
    print(paste0(param, ' mean Estimate: ', round(mean_est, 3)))
    print(paste0('95% CIs (mean): [', round(lower_ci_mean, 3), ' , ', round(upper_ci_mean, 3), ']'))
    print(paste0('mean bias: ', round(mean(bias, na.rm = TRUE), 3)))
    print(paste0('MAE: ', round(mean(MAE, na.rm = TRUE), 3)))
    print(paste0('mean sd: ', round(mean(sd, na.rm = TRUE), 3)))
    print(paste0('coverage: ', sum(coverage)))
    print(paste0('% coverage: ', coverage_pc))
    print(paste0(param, 'Effective Size (mean): ', round(mean_eff, 3)))
    
  }

#ADD GELMAN.DIAG
#gelman.diag(as.mcmc.list(c(1,1,2,1,2,2,1,3,3,1,1,1,1,1,1,1)))
#df3 = cbind(df1, df2)