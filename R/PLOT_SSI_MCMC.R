#PLOT SSI MCMC GRID
PLOT_SSI_MCMC <- function(epidemic_data, mcmc_output,
                               n_mcmc, cex = 1.8, RESULTS_FOLDER = '~/Github/computing/mcmc/SSI/',
                               PRIOR = TRUE, PDF = FALSE,
                               sim_vals = list(r0 = 2, k = 0.16), 
                               limits = list(r0 =  c(1.7, 2.3), k = c(0.075, 0.25)),
                               mcmc_specs = list(burn_in_pc = 0.2,
                                                 thinning_factor = 10)){
  
  #PLOT
  par(mfrow=c(4,2))
  par(mar = rep(4.5, 4), xpd = TRUE)

  #MODEL
  FLAGS_MODELS = GET_FLAGS_MODELS(SSI = TRUE)
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  MODEL_COLORS = GET_MODEL_COLORS(); MODEL_COLOR = MODEL_COLORS[3]
  
  #PARAMS
  FLAG_PARAM1 = GET_PARAM(r0 = TRUE); list_labels1 = GET_PARAM_LABEL(FLAG_PARAM1, model)
  FLAG_PARAM2 = GET_PARAM(k = TRUE); list_labels2 = GET_PARAM_LABEL(FLAG_PARAM2, model)
  
  #PRIORS 
  PRIORS_USED =  GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSI() 
  
  #MCMC
  r0_mcmc = mcmc_output$ssi_params_matrix[,1]; r0_mcmc = unlist(r0_mcmc); r0_mcmc = r0_mcmc[!is.na(r0_mcmc)]
  k_mcmc = mcmc_output$ssi_params_matrix[,2]; k_mcmc = unlist(k_mcmc); k_mcmc = k_mcmc[!is.na(k_mcmc)]
  
  if(PDF){
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_mcmc_', time_stamp, '.pdf') 
    create_folder(RESULTS_FOLDER)
    pdf(paste0(RESULTS_FOLDER, pdf_file), width = 14.0, height = 11.0)
    
    #MARGIN
    par(mfrow=c(4,2))
    par(mar = rep(5, 4), xpd = TRUE)
  }
  
  #LIMITS
  
  #LIMITS
  #r0_min =  min(sim_vals$r0, min(r0_mcmc, na.rm = TRUE));  r0_max =  max(sim_vals$r0, max(r0_mcmc, na.rm = TRUE))
  #k_min = min(sim_vals$k, min(k_mcmc, na.rm = TRUE)); k_max = max(sim_vals$k, max(k_mcmc, na.rm = TRUE))
  #r0_lim = c(r0_min, r0_max);  k_lim = c(k_min, k_max)
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #i. EPIDEMIC DATA
  PLOT_SIM_DATA(epidemic_data, FLAGS_MODELS)
  
  #ib. MARGINALS
  plot(r0_mcmc, k_mcmc,
       xlab = list_labels1$lab, ylab = list_labels2$lab,
       main = bquote(paste(italic(R[0]), ' vs ', .(list_labels2$lab))),
       col = MODEL_COLOR,
       cex.lab= cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
  
  #ii. TRACES 
  PLOT_MCMC_TRACE(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$r0, length(r0_mcmc), sim_vals$r0, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(k_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$k, length(k_mcmc), sim_vals$k, col = 'black', lwd = 2)
  
  #iii. HISTOGRAMS
  
  #HIST R0
  r0_mcmc_hist = subset(r0_mcmc, r0_mcmc > limits$r0[1])
  PLOT_MCMC_HIST(r0_mcmc_hist, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, xlim = limits$r0)
  PLOT_PRIOR_DIST(FLAG_PARAM1, r0_mcmc, limits$r0)
  segments(sim_vals$r0, 0, sim_vals$r0, 15, col = 'black', lwd = 2)

  #HIST K
  k_mcmc_hist = subset(k_mcmc, k_mcmc > limits$k[1])
  PLOT_MCMC_HIST(k_mcmc_hist, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, xlim = limits$k)
  PLOT_PRIOR_DIST(FLAG_PARAM2, k_mcmc, limits$k)
  segments(sim_vals$k, 0, sim_vals$k, 1.5, col = 'black', lwd = 2)

  #iv. CUMULATIVE MEANS 
  PLOT_CUMULATIVE_MEAN(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, ylim = limits$r0)
  segments(0, sim_vals$r0, length(r0_mcmc), sim_vals$r0, col = 'black', lwd = 2)
  
  PLOT_CUMULATIVE_MEAN(k_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, ylim = limits$k)
  segments(0, sim_vals$k, length(k_mcmc), sim_vals$k, col = 'black', lwd = 2)
  
  if(PDF){
    dev.off()
  }
  
  df_results <- data.frame(
    n_mcmc = n_mcmc,
    r0_sim = sim_vals$r0[1],
    r0_mean_mcmc = round(mean(r0_mcmc), 2), 
    r0_lo_95_cred_int = round(get_lower_ci(r0_mcmc), 2), 
    r0_up_95_cred_int = round(get_upper_ci(r0_mcmc), 2), 
    k_sim = sim_vals$k, 
    k_mean_mcmc = round(mean(k_mcmc), 2),
    k_lo_95_cred_int = round(get_lower_ci(k_mcmc), 2),
    k_up_95_cred_int = round(get_upper_ci(k_mcmc), 2), 
    accept_rate = round(mcmc_output$accept_rate, 2),
    r0_es = round(effectiveSize(as.mcmc(r0_mcmc))[[1]], 2),
    k_es = round(effectiveSize(as.mcmc(k_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) #format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
  
  #print(df_results)
  
  return(df_results)
}
