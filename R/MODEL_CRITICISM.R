#MODEL CRITICISM


GET_P_VALUE <- function(column_true_val, column_summary_stat) {
  'Get p values - comparing  summary stat columns to true value'
  
  #Final val
  num_iters = length(column_summary_stat)# - 1]
  #print(paste0('num_iters = ', num_iters))
  #P value
  prop_lt = length(which(column_summary_stat < column_true_val))/num_iters + 0.5*(length(which(column_summary_stat == column_true_val)))/num_iters
  #print(paste0('prop_lt = ', prop_lt))
  prop_gt = length(which(column_summary_stat > column_true_val))/num_iters + 0.5*(length(which(column_summary_stat == column_true_val)))/num_iters
  #print(paste0('prop_gt = ', prop_gt))
  pvalue = min(prop_lt, prop_gt)
  #print(paste0('pvalue 1 = ', pvalue))
  
  pvalue = 2*pvalue
  #print(paste0('pvalue 2 = ', pvalue))
  
  return(pvalue)
  
}


GET_SUMMARY_STATS <-function(epidemic_data) {
  
  mid_point = length(epidemic_data)/2; end_point = length(epidemic_data)
  
  summary_stats = c(
    sum_infects = sum(epidemic_data),
    sd_infect = sd(epidemic_data),
  sum_1st_half  = sum(epidemic_data[1:mid_point]),
  sum_2nd_half =  sum(epidemic_data[(mid_point+1):end_point]),
  
  median_infect = median(epidemic_data),
  max_infect = max(epidemic_data),
  
  #Differences
  max_dif = max((diff(epidemic_data))), #Change from absolute difference
  med_dif = median(diff(epidemic_data)),
  max_dif2nd = max(diff(diff(epidemic_data))),
  med_dif2nd =  median(diff(diff(epidemic_data))))

return(summary_stats )
}


#***************
#* GET_SIM_DATA

GET_SIM_DATA <- function(epidemic_data, mcmc_output, iter, FLAGS_MODELS){
  
  #PARAMS
  num_days = length(epidemic_data)
  
  if(FLAGS_MODELS$BASELINE){
    
    replicate_data = SIMULATE_EPI_BASELINE(num_days = num_days)
    
  } else if (FLAGS_MODELS$SSE){
    
    r0 = mcmc_output$sse_params_matrix[iter, 1]
    k = mcmc_output$sse_params_matrix[iter, 2]
    replicate_data = SIMULATE_EPI_SSE(num_days = num_days, r0 = r0, k = k)
    
  } else if (FLAGS_MODELS$SSI){
    
    r0 = mcmc_output$ssi_params_matrix[iter, 1]
    k = mcmc_output$ssi_params_matrix[iter, 2]
    epidemic_data_vec = SIMULATE_EPI_SSI(num_days = num_days, r0 = r0, k = k)
    replicate_data = epidemic_data_vec$epidemic_data
    
  } else if (FLAGS_MODELS$SSEB){
    
    r0 = mcmc_sseb$r0_vec[iter]
    alpha = mcmc_sseb$alpha_vec[iter]
    beta = mcmc_sseb$beta_vec[iter]
    replicate_data = SIMULATE_EPI_SSEB(num_days = num_days, r0 = r0, 
                                       alpha = alpha, beta = beta)
    
  } else if (FLAGS_MODELS$SSIB){
    
    r0 = mcmc_output$ssib_params_matrix[iter, 1]
    a = mcmc_output$ssib_params_matrix[iter, 2]
    b = mcmc_output$ssib_params_matrix[iter, 3]
    replicate_data = SIMULATE_EPI_SSIB(num_days = num_days, r0 = r0, a = a, b = b)
    
  } 
  
  return(replicate_data)
}

RUN_MODEL_CRITICISM <- function(epidemic_data, mcmc_output, 
                                FLAGS_MODELS = list(BASELINE = FALSE, 
                                                    SSE = FALSE, SSI = FALSE, SSEB = FALSE,
                                                    SSIB = FALSE), num_samples = 5000){
  
  #GET SUMMARY STATS
  true_summary_stats <- GET_SUMMARY_STATS(epidemic_data) 
  stat_names <- names(true_summary_stats)
  num_stats <- length(true_summary_stats)
  
  #REPLICATE SUMMARY STATS - MATRIX 
  replicate_summary_stats <- data.frame(matrix(ncol = num_stats, nrow = num_samples))
  colnames(replicate_summary_stats) <- stat_names
    
    for (iter in 1:num_samples){
      if(iter %% 500 == 0){
        print(iter)
      }
      replicate_data = GET_SIM_DATA(epidemic_data, mcmc_output, iter, FLAGS_MODELS)
      replicate_summary_stats[iter, ] <- GET_SUMMARY_STATS(replicate_data)
    }
  
  # Calculate p-values for each summary statistic
  # Assuming stat_names is a vector of summary statistic names
  p_values <- sapply(stat_names, function(stat) {
    GET_P_VALUE(true_summary_stats[stat], replicate_summary_stats[[stat]])
  })
  
  print(stat_names)
  print(p_values)
  # Convert the list to a data frame
  #p_values_df <- as.data.frame(p_values)
  #colnames(p_values_df) <- stat_names
  #p_values_matrix <- as.matrix(p_values_df)

  #print(p_values_matrix)
  
  return(p_values)
}



# Print the p-values

