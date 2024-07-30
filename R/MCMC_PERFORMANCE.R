#*********************************
#*
#* MCMC PERFORMANCE METRICS
#* 
#* **********************************
#' @export 
SIM_PERFORMANCE <- function(df_results, FLAG_PARAM = GET_PARAM(r0 = TRUE),
                            SSEB = FALSE, PRIOR = FALSE){
  
  #Param
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  print(paste0(param, ' Performance metrics'))
  num_runs = length(df_results$true_r0)
  
  #Bias, MAE, coverage 
  true_val = mean(unlist(df_results[paste0('true_', param)]))
  min_true = min(unlist(df_results[paste0('true_', param)]))
  max_true = max(unlist(df_results[paste0('true_', param)]))
  
  #Mean
  mean_est = mean(unlist(df_results[paste0('mean_', param)]))
  lower_ci_mean = mean(unlist(df_results[paste0('lower_ci_', param)]))
  upper_ci_mean = mean(unlist(df_results[paste0('upper_ci_', param)]))
  MAE = unlist(as.vector(abs(df_results[paste0('true_', param)] - df_results[paste0('mean_', param)])))
  bias = unlist(as.vector(df_results[paste0('mean_', param)] - df_results[paste0('true_', param)])) 
  sd = unlist(as.vector(df_results[paste0('sd_', param)]))
  coverage = unlist(as.vector(df_results[paste0('coverage_', param)]))
  coverage_pc = 100*(sum(coverage)/num_runs)
  
  #Results
  print(paste0(param, ' True value (mean): ', round(true_val, 3)))
  print(paste0(param, ' True value range: [', round(min_true, 3), ',', round(max_true, 3), ']'))
  print('***** Estimates (Mean): *****')
  print(paste0(param, ' Mean Estimate: ', round(mean_est, 3)))
  print(paste0('95% CIs (Mean): [', round(lower_ci_mean, 3), ' , ', round(upper_ci_mean, 3), ']'))
  print(paste0('Bias (mean): ', round(mean(bias, na.rm = TRUE), 5)))
  print(paste0('MAE: ', round(mean(MAE, na.rm = TRUE), 3)))
  print(paste0('sd (mean): ', round(mean(sd, na.rm = TRUE), 3)))
  print(paste0('Coverage: ', sum(coverage)))
  print(paste0('% Coverage (Mean): ', coverage_pc))
  
  #ACCEPT RATE
  if(SSEB){
    mean_accept_rate = mean(unlist(df_results[paste0('accept_rate_', param)]))
    print(paste0('accept rate ', param, ':', round(mean_accept_rate, 3))) 
  } else {
    accept_rate = mean(unlist((df_results['accept_rate'])))
    print(paste0('accept rate: ', round(accept_rate, 3))) 
  }
  
  if(PRIOR){
    print('no eff size')
  } else {
    mean_eff = mean(unlist(df_results[paste0('eff_size_', param)]))
    print(paste0(param, 'Effective Sample Size (mean): ', round(mean_eff, 3)))
  }
  
}


#SINGLE MCMC 
SIM_PERFORMANCE_MCMC <- function(mcmc_vec, true_val, FLAG_PARAM = GET_PARAM(r0 = TRUE),
                                 FLAGS_MODELS = GET_FLAGS_MODELS(BASELINE = TRUE)){
  
  #PARAM
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  #PERFORMANCE METRICS
  mean_est_posterior = round(mean(mcmc_vec), 3)
  bias_mcmc = round((mean_est_posterior - true_val), 3)
  lower_95_ci = round(get_lower_ci(mcmc_vec), 2)
  upper_95_ci = round(get_upper_ci(mcmc_vec), 2) 
  ess_mcmc = round(effectiveSize(as.mcmc(mcmc_vec))[[1]], 2)
  
  #RESULTS
  print(paste0(param, ' ', model, '. True value: ', true_val)) 
  print('***** Estimates: *****')
  print(paste0(param, ' Mean Estimate: ', mean_est_posterior))
  print(paste0('Bias: ', bias_mcmc))
  print(paste0('95% CIs (Mean): [', lower_95_ci, ' , ', upper_95_ci, ']'))
  print(paste0('Effective Sample Size: ', ess_mcmc))
  
}



