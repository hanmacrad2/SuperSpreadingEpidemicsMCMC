#PLOT SSI MCMC GRID
PLOT_SSIB_MCMC_PARAM_I <- function(epidemic_data, mcmc_output,
                           n_mcmc, list_epi_data, cex = 1.8, RESULTS_FOLDER = '~/Github/computing/mcmc/SSIB/SSIB_parameterisation_I',
                           PDF = FALSE, sim_vals = list(a = 0.5, b= 0.1, c = 10)){
  
  #PLOT
  par(mfrow=c(4,3))
  par(mar = rep(4.5, 4), xpd = TRUE)
  
  #MODEL
  FLAGS_MODELS = GET_FLAGS_MODELS(SSIB = TRUE)
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  MODEL_COLORS = GET_MODEL_COLORS(); MODEL_COLOR = MODEL_COLORS[5]
  
  #PARAMS
  FLAG_PARAM1 = GET_PARAM(a = TRUE); list_labels1 = GET_PARAM_LABEL_SSIB_I(FLAG_PARAM1, model, sim_vals$a)
  FLAG_PARAM2 = GET_PARAM(b = TRUE); list_labels2 = GET_PARAM_LABEL_SSIB_I(FLAG_PARAM2, model, sim_vals$b)
  FLAG_PARAM3 = GET_PARAM(c = TRUE); list_labels3 = GET_PARAM_LABEL_SSIB_I(FLAG_PARAM3, model, sim_vals$c)
  
  #MCMC
  a_mcmc = mcmc_output$a_vec; b_mcmc = mcmc_output$b_vec; c_mcmc = mcmc_output$c_vec 
  
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
  a_lim = c(min(min(a_mcmc), sim_vals$a), max(max(a_mcmc), sim_vals$a))  
  b_lim = c(min(min(b_mcmc), sim_vals$b), max(max(b_mcmc), sim_vals$b))   
  c_lim = c(min(min(c_mcmc), sim_vals$c), max(max(c_mcmc), sim_vals$c)) #c(min(c_mcmc), max(c_mcmc))
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #i. EPIDEMIC DATA
  PLOT_SIM_DATA_SSIB_I(epidemic_data, sim_vals)
  PLOT_LOG_LIKELIHOOD(mcmc_output$log_like_vec, FLAGS_MODELS, n_mcmc) # cex = cex)
  PLOT_SSIB_DATA(mcmc_output, list_epi_data, cex)
  #plot.ts(0, xlab = '', ylab = '',  xaxt = "n", yaxt = "n")
  
  #ii. TRACES 
  PLOT_MCMC_TRACE(a_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$a, length(a_mcmc), sim_vals$a, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(b_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$b, length(b_mcmc), sim_vals$b, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(c_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$c, length(c_mcmc), sim_vals$c, col = 'black', lwd = 2)
  
  #iii. HISTOGRAMS
  
  #HIST a
  a_mcmc_hist = subset(a_mcmc, a_mcmc > a_lim[1])
  PLOT_MCMC_HIST_SSIB_I(a_mcmc_hist, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, xlim = a_lim)
  PLOT_PRIOR_DIST_SSIB_I(FLAG_PARAM1, a_mcmc, a_lim)
  segments(sim_vals$a, 0, sim_vals$a, 12, col = 'black', lwd = 2)
  
  #HIST b
  PLOT_MCMC_HIST_SSIB_I(b_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, xlim = b_lim)
  PLOT_PRIOR_DIST_SSIB_I(FLAG_PARAM2, b_mcmc, b_lim)
  segments(sim_vals$b, 0, sim_vals$b, 100, col = 'black', lwd = 2)
  
  #HIST c
  #c_mcmc_hist = subset(c_mcmc, c_mcmc > 5)
  PLOT_MCMC_HIST_SSIB_I(c_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex, xlim = c_lim)
  PLOT_PRIOR_DIST_SSIB_I(FLAG_PARAM3, c_mcmc, c_lim)
  segments(sim_vals$c, 0, sim_vals$c, 0.7, col = 'black', lwd = 2)
  
  #iv. CUMULATIVE MEANS 
  PLOT_CUMULATIVE_MEAN(a_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, ylim = a_lim)
  segments(0, sim_vals$a, length(a_mcmc), sim_vals$a, col = 'black', lwd = 2)
  
  PLOT_CUMULATIVE_MEAN(b_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, ylim = b_lim)
  segments(0, sim_vals$b, length(b_mcmc), sim_vals$b, col = 'black', lwd = 2, ylim = b_lim)
  
  PLOT_CUMULATIVE_MEAN(c_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex, ylim = c_lim)
  segments(0, sim_vals$c, length(c_mcmc), sim_vals$c, col = 'black', lwd = 2)
  
  if(PDF){
    dev.off()
  }
  
  df_results <- data.frame(
    n_mcmc = n_mcmc,
    a_sim = sim_vals$a,
    a_mean_mcmc = round(mean(a_mcmc), 2), 
    a_lo_95_cred_int = round(get_lower_ci(a_mcmc), 2), 
    a_up_95_cred_int = round(get_upper_ci(a_mcmc), 2), 
    b_sim = sim_vals$b, 
    b_mean_mcmc = round(mean(b_mcmc), 2),
    b_lo_95_cred_int = round(get_lower_ci(b_mcmc), 2),
    b_up_95_cred_int = round(get_upper_ci(b_mcmc), 2),
    c_sim = sim_vals$c, 
    c_mean_mcmc = round(mean(c_mcmc), 2),
    c_lo_95_cred_int = round(get_lower_ci(c_mcmc), 2),
    c_up_95_cred_int = round(get_upper_ci(c_mcmc), 2), 
    #accept_rate_r0 = round(mcmc_output$list_accept_rates$accept_rate2, 2),
    #accept_rate_a = round(mcmc_output$list_accept_rates$accept_rate1, 2),
    #saccept_rate_b = round(mcmc_output$list_accept_rates$accept_rate3, 2),
    a_es = round(effectiveSize(as.mcmc(a_mcmc))[[1]], 2),
    b_es = round(effectiveSize(as.mcmc(b_mcmc))[[1]], 2),
    c_es = round(effectiveSize(as.mcmc(c_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) 
  
  #print(df_results)
  
  return(df_results)
}


#PLOT SIM DATA SSIB I
PLOT_SIM_DATA_SSIB_I <- function(epidemic_data, sim_vals, main_font = 2.7, cex = 1.6){
  
  data_title = bquote(paste("SSIB simulated data: a = ", .(sim_vals$a),
  " b = ", .(sim_vals$b), " c = ", .(sim_vals$c)))
  #PLOT
  plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
       xlab = 'Time', ylab = 'Daily Infection count',
       main = data_title,
       lwd = 2.5, #3.5,
       cex.lab= cex, cex.axis= cex, cex.main= cex, cex.sub=cex)
  
}


#PARAM LABELS FOR PLOT
GET_PARAM_LABEL_SSIB_I <- function(FLAG_PARAM, model, sim_val) { #model or FLAG_MODELS
  
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  #list labels
  list_labels = list(lab = bquote(paste(.(param))), 
                     main_trace =  bquote(paste(.(param), " Trace - ", .(model)~ "model")),
                     main_hist =  bquote(paste(.(param), " Posterior - ", .(model)~ "model")),
                     main_mean2 =  bquote(paste(.(param), " Cumulative mean - ", .(model)~ "model")),
                     main_hist_prior = bquote(paste(.(param), " Posterior. Prior: Exponential(1)")),
                     main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated value = ", .(sim_val))))
  
  if(FLAG_PARAM$c){
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: 1 + Gamma(3,3)"))
  }
  
  return(list_labels)
  
}
