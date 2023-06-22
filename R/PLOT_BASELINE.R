#PLOT R0 RESELTS
#' Grid Plot of Baseline MCMC results for R0
#'
#'Grid Plot of MCMC Results for Super-Spreading Epidemic Models
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_output mcmc samples from mcmc sampler/algorithm
#' 
PLOT_BASELINE_R0_MCMC <- function(epidemic_data, mcmc_output, true_r0 = 1.6,
                                  true_loglike = 0, data_type = 'Baseline', ADAPTIVE = TRUE) { #sim_data, mcmc_output, true_r0, time_elap, seed_count, model_type){
  
  #PLOT
  plot.new(); par(mfrow=c(2,3))
  
  #MCMC OUTPUT
  r0_mcmc = mcmc_output$r0_vec; r0_mcmc = unlist(r0_mcmc)
  #Limits
  ll_lim_min = min(true_loglike, min(mcmc_output$log_like_vec, na.rm = TRUE))
  ll_lim_max = max(true_loglike, max(mcmc_output$log_like_vec, na.rm = TRUE))
  
  #***********
  #* PLOTS *
  
  #i. EPIDEMIC DATA
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count', 
          main = paste0(data_type, " Epidemic data"), 
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MEAN PLOTS
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  plot(seq_along(r0_mean), r0_mean,
               xlab = 'Time', ylab = 'R0', main = paste('R0 MCMC Mean'),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #iii. RO SAMPLES - HISTOGRAM
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', 
       main = paste('R0 total MCMC samples'), 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  #iv. MCMC TRACE PLOTS
  plot.ts(r0_mcmc, ylab = 'R0', 
          main = paste('R0 MCMC Baseline Model'),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #PLOT LOG-LIKE
  plot.ts(mcmc_output$log_like_vec, ylab = 'Log-likelihood', 
          main = paste('Log-likelihood MCMC'),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
          ylim = c(ll_lim_min, ll_lim_max))
  abline(h = true_loglike, col = 'orange', lwd = 2)

  plot.ts(mcmc_output$sigma_vec, ylab = 'R0', 
          main = paste('Sigma Adaptive MCMC'),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #FINAL MEAN STATS 
  r0_mcmc_mean = round(mean(r0_mcmc[length(r0_mcmc)/2:length(r0_mcmc)]), 2)
  
  #Results
  df_results <- data.frame(
    R0_mean_mc = r0_mcmc_mean,
    accept_rate_r0 = round(mcmc_output[[4]],2)) #time_elap = round(time_elap,2)) 
  
  return(df_results)
  
}