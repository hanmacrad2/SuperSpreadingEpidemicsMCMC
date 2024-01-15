#PLOT R0 RESELTS
#' Grid Plot of Baseline MCMC results for R0
#'
#'Grid Plot of MCMC Results for Super-Spreading Epidemic Models
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_output mcmc samples from mcmc sampler/algorithm
#' 
PLOT_BASELINE_MCMC <- function(epidemic_data, mcmc_output,
                               r0_sim = 2.0, true_loglike = 0, cex = 1.6,
                                  PDF = TRUE, ADAPTIVE = TRUE,
                                  PRIORS = list(EXP = FALSE, UNIF = FALSE, GAMMA = FALSE)) { #sim_data, mcmc_output, r0_sim, time_elap, seed_count, model_type){
  
  #PLOT
  #plot.new()
  #MODEL
  FLAGS_MODELS = GET_FLAGS_MODELS(BASELINE = TRUE)
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  MODEL_COLORS = GET_MODEL_COLORS(); MODEL_COLOR = MODEL_COLORS[1]
  
  if(PDF){
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_mcmc_', time_stamp, '.pdf') 
    RESULTS_FOLDER = '~/Github/Results/MCMC/Baseline/'
    pdf(paste0(RESULTS_FOLDER, pdf_file), width = 12.0, height = 8.0)
    
    #MARGIN
    par(mfrow=c(3,2))
    par(mar = rep(5, 4), xpd = TRUE)
    #par(mar=c(5.2, 4.8, 3.0, 19.45), xpd=TRUE) #Margins; bottom, left, top, right
    #par(mar= rep(5.0, 4), xpd=TRUE) 
  }
  
  #MCMC OUTPUT
  r0_mcmc = mcmc_output$r0_vec; r0_mcmc = unlist(r0_mcmc)
  #Limits
  ll_lim_min = min(true_loglike, min(mcmc_output$log_like_vec, na.rm = TRUE))
  ll_lim_max = max(true_loglike, max(mcmc_output$log_like_vec, na.rm = TRUE))
  
  #***********
  #* PLOTS *
  
  #i. EPIDEMIC DATAs
  PLOT_SIM_DATA(epidemic_data, FLAGS_MODELS)
  
  #i. MCMC TRACE PLOTS
  PLOT_R0_TRACE(r0_mcmc, FLAGS_MODELS, MODEL_COLOR, cex = cex)
  segments(0, r0_sim, length(r0_mcmc), r0_sim, col = 'black', lwd = 2)
  #abline(h = r0_sim, col = 'black', lwd = 2, xlim = c(0, max(r0_mcmc)))
  
  #iv. Cumulative Mean
  PLOT_CUMULATIVE_MEAN(r0_mcmc, FLAGS_MODELS, MODEL_COLOR, cex = cex)
  segments(0, r0_sim, length(r0_mcmc), r0_sim, col = 'black', lwd = 2)
  
  #iii. RO SAMPLES - HISTOGRAM
  PLOT_HISTOGRAM(r0_mcmc, FLAGS_MODELS, MODEL_COLOR, cex = cex)
  segments(r0_sim, 0, r0_sim, 10, col = 'black', lwd = 2)
  #lines(c(0, 9), ylim = c(0, 9), col = 'black', lwd = 2)
  #abline(v = r0_sim, col = 'black', lwd = 2, ylim = c(0, 9))
  
  #v. PLOT LOG-LIKE
  PLOT_LOG_LIKELIHOOD(mcmc_output$log_like_vec, FLAGS_MODELS)
  
  #vi. PLOT SIGMA
  PLOT_SIGMA(mcmc_output$sigma_vec)
  
  #PRIOR TITLES
  if(PRIORS$EXP){
    prior_title = 'Exponential(1)' 
  } else if (PRIORS$UNIF){
    prior_title = 'Uniform(0,10)' 
  } else if (PRIORS$GAMMA){
    prior_title = 'Gamma(1,5)' 
  }
  
  #PRIORS
  if(PRIORS$EXP){
    mr0_prior = 'exp(1)' #paste0('exp(', list_priors$r0[1], ')')
    xseq_r0 = seq(0, 3, length.out = 500)
    dr0e = dexp(xseq_r0, 1) #list_priors$r0[1])
    lines(xseq_r0, dr0e, type = 'l', lwd = 2)
    
  } else if(PRIORS$UNIF){
    mr0_prior = 'exp(1)' #paste0('exp(', list_priors$r0[1], ')')
    xseq_r0 = seq(0, 15, length.out = 1000)
    dr0e = dunif(xseq_r0, min = 0, max = 10) #list_priors$r0[1])
    lines(xseq_r0, dr0e, type = 'l', lwd = 2)
  }
  
  #Results
  df_results <- data.frame(
    R0_mean_mc = round(mean(r0_mcmc[length(r0_mcmc)/2:length(r0_mcmc)]), 2),
    accept_rate_r0 = round(mcmc_output[[4]],2)) #time_elap = round(time_elap,2)) 
  
  if(PDF){
    dev.off()
  }
  return(df_results)
  
}