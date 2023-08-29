#SSE SIMULATIONS
par(mfrow = c(3,2))
PLOT_BOUNDS <-function(num_days = 50, R0 = 1.5, k = 0.2){
  
  #MATRIX OF SIMULATIONS
  list_matrix_sims = GET_MATRIX_SIM_SSE(R0 = R0, k = k)
  matrix_sim_data = list_matrix_sims$matrix_sim_temp
  model_detail = list_matrix_sims$model_detail
  
  #BOUNDS 
  mean_est <- apply(matrix_sim_data, 2, mean) |> unlist()
  upper_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.975) |> unlist()
  lower_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.025) |> unlist()
  posterior_zig_zag = apply(matrix_sim_data, 1, zigzag)
  
  #PLOT
  ylim = c(min(lower_bounds), max(upper_bounds))
  plot(1:num_days, mean_est, type = 'l', ylim = ylim,
       main = model_detail, xlab = 'time', ylab = 'Infection count')
  lines(1:num_days, upper_bounds, col = 'red', lwd = 2) 
  lines(1:num_days, lower_bounds, col = 'red', lwd = 2) 
  
  #HISTOGRAM OF ZIG-ZAG 
  hist(posterior_zig_zag, breaks = 50, 
       main = paste0('Zig-zag of sim data: sum(abs(diff(data))). ', model_detail), xlab = 'zig zag of sim data')
  abline(v = quantile(posterior_zig_zag, probs = c(0.025, 0.975)), col = 'red', lwd = 2)
  abline(v = mean(posterior_zig_zag), col = 'black', lwd = 2)
}


GET_MATRIX_SIM_SSE <- function(n_sample_repeats = 1000, num_days = 50,
                               model_type = 'SSE',
                               R0 = 1.5, k = 0.2) {
  #PARAMS
  matrix_sim_temp = matrix(nrow = n_sample_repeats, ncol = num_days)
  model_detail = paste0(model_type, ' model. Mean(black) 95% CIs(red) of ', n_sample_repeats, ' sims. R0: ', R0, ' k: ', k)
  print(paste0('R0:, ', R0)); print(paste0('k:, ', k))
  
  #SIMULATE
  for (i in 1:n_sample_repeats){
    
    sim_data_sse = SIMULATE_EPI_SSE(R0 = R0, k = k)
    matrix_sim_temp[i, ] = sim_data_sse
    
  }
 
  return(list(matrix_sim_temp = matrix_sim_temp, model_detail = model_detail)) 
}


#APPLY
#PLOT_BOUNDS(k = 10)
