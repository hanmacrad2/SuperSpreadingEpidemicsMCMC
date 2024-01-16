#PLOT SSI MCMC GRID
PLOT_SSEB_MCMC <- function(epidemic_data, mcmc_output,
                               n_mcmc, cex = 1.8, RESULTS_FOLDER = '~/Github/computing/mcmc/SSEB/',
                               PRIOR = FALSE, PDF = FALSE,
                               sim_vals = list(r0 = 2, alpha= 0.5, beta = 10), 
                               mcmc_specs = list(burn_in_pc = 0.2,
                                                 thinning_factor = 10)){
  
  #PLOT
  par(mfrow=c(4,3))
  par(mar = rep(4.5, 4), xpd = TRUE)
  
  #MODEL
  FLAGS_MODELS = GET_FLAGS_MODELS(SSEB = TRUE)
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  MODEL_COLORS = GET_MODEL_COLORS(); MODEL_COLOR = MODEL_COLORS[4]
  
  #PARAMS
  FLAG_PARAM1 = GET_PARAM(r0 = TRUE); list_labels1 = GET_PARAM_LABEL(FLAG_PARAM1, model)
  FLAG_PARAM2 = GET_PARAM(alpha = TRUE); list_labels2 = GET_PARAM_LABEL(FLAG_PARAM2, model)
  FLAG_PARAM3 = GET_PARAM(beta = TRUE); list_labels3 = GET_PARAM_LABEL(FLAG_PARAM3, model)
  
  #PRIORS 
  PRIORS_USED =  GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSEB() 
  
  #MCMC
  r0_mcmc = mcmc_output$r0_vec; alpha_mcmc = mcmc_output$alpha_vec; beta_mcmc = mcmc_output$beta_vec 
  
  if(PDF){
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_mcmc_', time_stamp, '.pdf') 
    create_folder(RESULTS_FOLDER)
    pdf(paste0(RESULTS_FOLDER, pdf_file), width = 14.0, height = 11.0)
    
    #MARGIN
    par(mfrow=c(4,3))
    par(mar = rep(5, 4), xpd = TRUE)
  }
  
  #LIMITS
  r0_lim = c(0.9, 3.0); alpha_lim = c(0,1); beta_lim = c(5, 15)
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #i. EPIDEMIC DATA
  plot.ts(0, xlab = '', ylab = '',  xaxt = "n", yaxt = "n")
  
  PLOT_SIM_DATA(epidemic_data, FLAGS_MODELS)
  plot.ts(0, xlab = '', ylab = '',  xaxt = "n", yaxt = "n")

  #ii. TRACES 
  PLOT_MCMC_TRACE(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$r0, length(r0_mcmc), sim_vals$r0, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(alpha_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$alpha, length(alpha_mcmc), sim_vals$alpha, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(beta_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$beta, length(beta_mcmc), sim_vals$beta, col = 'black', lwd = 2)
  
  #iii. HISTOGRAMS
  
  #HIST R0
  r0_mcmc_hist = subset(r0_mcmc, r0_mcmc > 0.95)
  PLOT_MCMC_HIST(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, xlim = r0_lim)
  PLOT_PRIOR_DIST(FLAG_PARAM1, r0_mcmc, r0_lim)
  segments(sim_vals$r0, 0, sim_vals$r0, 1.0, col = 'black', lwd = 2)
  
  #HIST ALPHA
  PLOT_MCMC_HIST(alpha_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, xlim = alpha_lim)
  PLOT_PRIOR_DIST(FLAG_PARAM2, alpha_mcmc, alpha_lim)
  segments(sim_vals$alpha, 0, sim_vals$alpha, 4, col = 'black', lwd = 2)

  #HIST BETA
  beta_mcmc_hist = subset(beta_mcmc, beta_mcmc > 5)
  PLOT_MCMC_HIST(beta_mcmc_hist, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex, xlim = beta_lim)
  PLOT_PRIOR_DIST(FLAG_PARAM3, beta_mcmc, beta_lim)
  segments(sim_vals$beta, 0, sim_vals$beta, 0.25, col = 'black', lwd = 2)

  #iv. CUMULATIVE MEANS 
  PLOT_CUMULATIVE_MEAN(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, ylim = r0_lim)
  segments(0, sim_vals$r0, length(r0_mcmc), sim_vals$r0, col = 'black', lwd = 2)
  
  PLOT_CUMULATIVE_MEAN(alpha_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, ylim = alpha_lim)
  segments(0, sim_vals$alpha, length(alpha_mcmc), sim_vals$alpha, col = 'black', lwd = 2, ylim = alpha_lim)
  
  PLOT_CUMULATIVE_MEAN(beta_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex, ylim = beta_lim)
  segments(0, sim_vals$beta, length(beta_mcmc), sim_vals$beta, col = 'black', lwd = 2)
  
  if(PDF){
    dev.off()
  }
  
  df_results <- data.frame(
    n_mcmc = n_mcmc,
    r0_sim = sim_vals$r0[1],
    r0_mean_mcmc = round(mean(r0_mcmc), 2), 
    r0_lo_95_cred_int = round(get_lower_ci(r0_mcmc), 2), 
    r0_up_95_cred_int = round(get_upper_ci(r0_mcmc), 2), 
    alpha_sim = sim_vals$alpha, 
    alpha_mean_mcmc = round(mean(alpha_mcmc), 2),
    alpha_lo_95_cred_int = round(get_lower_ci(alpha_mcmc), 2),
    alpha_up_95_cred_int = round(get_upper_ci(alpha_mcmc), 2),
    beta_sim = sim_vals$beta, 
    beta_mean_mcmc = round(mean(beta_mcmc), 2),
    beta_lo_95_cred_int = round(get_lower_ci(beta_mcmc), 2),
    beta_up_95_cred_int = round(get_upper_ci(beta_mcmc), 2), 
    accept_rate_r0 = round(mcmc_output$list_accept_rates$accept_rate_r0, 2),
    accept_rate_alpha = round(mcmc_output$list_accept_rates$accept_rate_alpha, 2),
    accept_rate_beta = round(mcmc_output$list_accept_rates$accept_rate_beta, 2),
    r0_es = round(effectiveSize(as.mcmc(r0_mcmc))[[1]], 2),
    alpha_es = round(effectiveSize(as.mcmc(alpha_mcmc))[[1]], 2),
    beta_es = round(effectiveSize(as.mcmc(beta_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) 
  
  #print(df_results)
  
  return(df_results)
}
