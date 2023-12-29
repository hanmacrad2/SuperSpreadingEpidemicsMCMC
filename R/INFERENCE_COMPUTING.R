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
  
  return(row_k)
}

#PARAMTER INFERENCE
GET_PARAM_INFERENCE <- function(true_val, param_vec,
                                FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE,
                                                  alpha = FALSE, a = FALSE, 
                                                  beta = FALSE, b = FALSE)){
  
  # Get names of true parameters
  param_name <- names(FLAG_PARAM)[FLAG_PARAM]
  
  # Create a function to generate column names
  generate_column_names <- function(prefix, param_name) {
    
    return(paste(prefix, param_name, sep = "_"))
  }
  
  # Create column names for different statistics
  mean_prefix <- "mean"
  true_prefix <- "true"
  col_names <- c(
    generate_column_names("true", param_name),
    generate_column_names("mean", param_name),
    generate_column_names("sd", param_name),
    generate_column_names("lower_ci_", param_name),
    generate_column_names("upper_ci_", param_name),
    generate_column_names("coverage_", param_name),
    generate_column_names("esize_", param_name)
  )
  
  # Create a data frame with new column names
  row_result <- data.frame(
    !!col_names[1] := true_val,
    !!col_names[2] := mean(param_vec),
    !!col_names[3] := sd(param_vec),
    !!col_names[4] := get_lower_ci(param_vec),
    !!col_names[5] := get_upper_ci(param_vec),
    !!col_names[6] := GET_COVERAGE(param_vec),
    !!col_names[7] := effectiveSize(param_vec)
  )
  
  colnames(row_result) <- col_names
  
  return(row_result)
}

GET_FLAG_PARAM <- function(){
  
  FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE,
                    alpha = FALSE, a = FALSE, 
                    beta = FALSE, b = FALSE)
  
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
#* *************************************************************

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

#*******************************************
# 2. SSE MODEL
#*******************************************
#' @export 
INFER_SSE <- function(r0_val, k_val, PRIORS_USED, num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  epidemic_data = SIMULATE_EPI_SSE(num_days = num_days, r0 = r0_val, k = k_val)
  
  #MCMC
  mcmc_output = MCMC_INFER_SSE(epidemic_data, n_mcmc, PRIORS_USED)
  r0_vec = mcmc_output$sse_params_matrix[,1]
  k_vec = mcmc_output$sse_params_matrix[,2]
  
  #Row result
  r0_row = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  k_row = GET_INFER_K_ROW(k_val, k_vec)
  
  result_row <- cbind(r0_row, k_row)
  
  return(result_row)
}

#*******************************************
# 3. SSI MODEL
#*******************************************
#' @export 
INFER_SSI <- function(r0_val, k_val, PRIORS_USED, num_days = 50, n_mcmc = 40000){ #{40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  data_ssi = SIMULATE_EPI_SSI(num_days = num_days, r0 = r0_val, k = k_val)
  epidemic_data = data_ssi$epidemic_data
  
  #MCMC
  mcmc_output = MCMC_INFER_SSI(epidemic_data, n_mcmc, PRIORS_USED)
  r0_vec = mcmc_output$ssi_params_matrix[,1]
  k_vec = mcmc_output$ssi_params_matrix[,2]

  #Row result
  r0_row = GET_INFER_R0_ROW(r0_val, r0_vec, mcmc_output, epidemic_data)
  k_row = GET_INFER_K_ROW(k_val, k_vec)
  
  result_row <- cbind(r0_row, k_row)
  
  return(result_row)
}


#**********************************************
#* SSEB MODEL
#* ********************************************
INFER_SSEB <- function(r0_val, alpha_val, beta_val, PRIORS_USED,
                       num_days = 50, n_mcmc = 40000) {
  
  'Inference of baseline simulate data'
  cat(r0_val)
  n_mcmc = 100
  epidemic_data = SIMULATE_EPI_SSEB(num_days = num_days, r0 = r0_val,
                                alpha = alpha_val, beta = beta_val)
  
  #MCMC
  mcmc_output = MCMC_INFER_SSEB(epidemic_data, n_mcmc, PRIORS_USED)
  
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
  
  result_row = cbind(row_r0, row_alpha, row_beta)
  
  return(result_row)
}


#
#ADD GELMAN.DIAG
#gelman.diag(as.mcmc.list(c(1,1,2,1,2,2,1,3,3,1,1,1,1,1,1,1)))
#df3 = cbind(df1, df2)