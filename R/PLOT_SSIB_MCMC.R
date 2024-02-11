#PLOT SSI MCMC GRID
PLOT_SSIB_MCMC <- function(epidemic_data, mcmc_output,
                           n_mcmc, list_epi_data, cex = 1.8, RESULTS_FOLDER = '~/Github/computing/mcmc/SSIB/',
                           PDF = TRUE, JOINT = TRUE, PARAM_I = FALSE,
                           sim_vals = list(r0 = 2, a= 0.85, b = 10), 
                           mcmc_specs = list(burn_in_pc = 0.2,
                                             thinning_factor = 10)){
  
  #PLOT
  par(mfrow=c(4,3))
  par(mar = rep(4.5, 4), xpd = TRUE)
  
  #MODEL
  FLAGS_MODELS = GET_FLAGS_MODELS(SSIB = TRUE)
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  MODEL_COLORS = GET_MODEL_COLORS(); MODEL_COLOR = MODEL_COLORS[5]
  
  #PARAMS
  FLAG_PARAM1 = GET_PARAM(r0 = TRUE); list_labels1 = GET_PARAM_LABEL(FLAG_PARAM1, model)
  FLAG_PARAM2 = GET_PARAM(a = TRUE); list_labels2 = GET_PARAM_LABEL(FLAG_PARAM2, model)
  FLAG_PARAM3 = GET_PARAM(b = TRUE); list_labels3 = GET_PARAM_LABEL(FLAG_PARAM3, model)
  
  #MCMC
  r0_mcmc = mcmc_output$r0_vec; a_mcmc = mcmc_output$alpha_vec; b_mcmc = mcmc_output$b_vec 
  
  #MCMC
  if(JOINT){
    r0_mcmc = mcmc_output$ssib_params_matrix[,1]; r0_mcmc = unlist(r0_mcmc); r0_mcmc = r0_mcmc[!is.na(r0_mcmc)]
    a_mcmc = mcmc_output$ssib_params_matrix[,2]; a_mcmc = unlist(a_mcmc); a_mcmc = a_mcmc[!is.na(a_mcmc)]
    b_mcmc = mcmc_output$ssib_params_matrix[,3]; b_mcmc = unlist(b_mcmc); b_mcmc = b_mcmc[!is.na(b_mcmc)]
    
  } else if (PARAM_I){
    r0_mcmc = mcmc_output$a_vec; a_mcmc = mcmc_output$b_vec; b_mcmc = mcmc_output$c_vec  
    
    #PARAMS
    FLAG_PARAM1 = GET_PARAM(a = TRUE); list_labels1 = GET_PARAM_LABEL(FLAG_PARAM1, model)
    FLAG_PARAM2 = GET_PARAM(b = TRUE); list_labels2 = GET_PARAM_LABEL(FLAG_PARAM2, model)
    FLAG_PARAM3 = GET_PARAM(c = TRUE); list_labels3 = GET_PARAM_LABEL(FLAG_PARAM3, model)
  }
  
  if(PDF){
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_mcmc_', time_stamp, '.pdf') 
    RESULTS_FOLDER = paste0(RESULTS_FOLDER, 'plots/')
    create_folder(RESULTS_FOLDER)
    pdf(paste0(RESULTS_FOLDER, pdf_file), width = 14.0, height = 11.0)
    
    #MARGIN
    par(mfrow=c(4,3))
    par(mar = rep(5, 4), xpd = TRUE)
  }
  
  #LIMITS
  r0_lim = c(min(r0_mcmc), max(r0_mcmc)) #r0_lim = c(1.85, 2.2); 
  a_lim = c(0,1)
  b_lim = c(min(b_mcmc), max(b_mcmc))
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #i. EPIDEMIC DATA
  PLOT_SIM_DATA(epidemic_data, FLAGS_MODELS)
  PLOT_LOG_LIKELIHOOD(mcmc_output$log_like_vec, FLAGS_MODELS, n_mcmc) # cex = cex)
  PLOT_SSIB_DATA(mcmc_output, list_epi_data, cex)
  #plot.ts(0, xlab = '', ylab = '',  xaxt = "n", yaxt = "n")
  
  #ii. TRACES 
  PLOT_MCMC_TRACE(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$r0, length(r0_mcmc), sim_vals$r0, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(a_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$a, length(a_mcmc), sim_vals$a, col = 'black', lwd = 2)
  
  PLOT_MCMC_TRACE(b_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex)
  segments(0, sim_vals$b, length(b_mcmc), sim_vals$b, col = 'black', lwd = 2)
  
  #iii. HISTOGRAMS
  
  #HIST R0
  r0_mcmc_hist = subset(r0_mcmc, r0_mcmc > r0_lim[1])
  PLOT_MCMC_HIST(r0_mcmc_hist, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex) #, xlim = r0_lim)
  PLOT_PRIOR_DIST(FLAG_PARAM1, r0_mcmc, r0_lim)
  segments(sim_vals$r0, 0, sim_vals$r0, 5.0, col = 'black', lwd = 2)
  
  #HIST a
  PLOT_MCMC_HIST(a_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex) #, xlim = a_lim)
  PLOT_PRIOR_DIST(FLAG_PARAM2, a_mcmc, a_lim)
  segments(sim_vals$a, 0, sim_vals$a, 2, col = 'black', lwd = 2)
  
  #HIST b
  #b_mcmc_hist = subset(b_mcmc, b_mcmc > 5)
  PLOT_MCMC_HIST(b_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex) #, xlim = b_lim)
  PLOT_PRIOR_DIST(FLAG_PARAM3, b_mcmc, b_lim)
  segments(sim_vals$b, 0, sim_vals$b, 0.12, col = 'black', lwd = 2)
  
  #iv. CUMULATIVE MEANS 
  PLOT_CUMULATIVE_MEAN(r0_mcmc, FLAGS_MODELS, FLAG_PARAM1, MODEL_COLOR, cex = cex, ylim = r0_lim)
  segments(0, sim_vals$r0, length(r0_mcmc), sim_vals$r0, col = 'black', lwd = 2)
  
  PLOT_CUMULATIVE_MEAN(a_mcmc, FLAGS_MODELS, FLAG_PARAM2, MODEL_COLOR, cex = cex, ylim = a_lim)
  segments(0, sim_vals$a, length(a_mcmc), sim_vals$a, col = 'black', lwd = 2, ylim = a_lim)
  
  PLOT_CUMULATIVE_MEAN(b_mcmc, FLAGS_MODELS, FLAG_PARAM3, MODEL_COLOR, cex = cex, ylim = b_lim)
  segments(0, sim_vals$b, length(b_mcmc), sim_vals$b, col = 'black', lwd = 2)
  
  if(PDF){
    dev.off()
  }
  
  df_results <- data.frame(
    n_mcmc = n_mcmc,
    r0_sim = sim_vals$r0[1],
    r0_mean_mcmc = round(mean(r0_mcmc), 2), 
    r0_lo_95_cred_int = round(get_lower_ci(r0_mcmc), 2), 
    r0_up_95_cred_int = round(get_upper_ci(r0_mcmc), 2), 
    a_sim = sim_vals$a, 
    a_mean_mcmc = round(mean(a_mcmc), 2),
    a_lo_95_cred_int = round(get_lower_ci(a_mcmc), 2),
    a_up_95_cred_int = round(get_upper_ci(a_mcmc), 2),
    b_sim = sim_vals$b, 
    b_mean_mcmc = round(mean(b_mcmc), 2),
    b_lo_95_cred_int = round(get_lower_ci(b_mcmc), 2),
    b_up_95_cred_int = round(get_upper_ci(b_mcmc), 2), 
    #accept_rate_r0 = round(mcmc_output$list_accept_rates$accept_rate2, 2),
    #accept_rate_a = round(mcmc_output$list_accept_rates$accept_rate1, 2),
    #saccept_rate_b = round(mcmc_output$list_accept_rates$accept_rate3, 2),
    r0_es = round(effectiveSize(as.mcmc(r0_mcmc))[[1]], 2),
    a_es = round(effectiveSize(as.mcmc(a_mcmc))[[1]], 2),
    b_es = round(effectiveSize(as.mcmc(b_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) 
  
  #print(df_results)
  
  return(df_results)
}


#PLOT
PLOT_SSIB_DATA <- function(mcmc_ssib, list_data_ssib, cex, lwd = 2.0, 
                           col_ns = 'aquamarine', col_ss = 'orange'){
  
  #Plot TRUE
  non_ss = list_data_ssib$non_ss
  ss = list_data_ssib$ss
  ns_inf = mcmc_ssib$data[[1]]
  ss_inf = mcmc_ssib$data[[2]]
  title = bquote('True vs Inferred: Non SSIs (aqua), SSI (or)')
  
  #nsI TRUE
  plot(seq_along(non_ss), non_ss,  type = 'l',
       xlab = 'Time', ylab = 'Daily infection count',
       main = title,
       col = col_ns,
       lwd = lwd, 
       cex.lab= cex, cex.axis= cex, cex.main= cex, cex.sub=cex)
  
  #nsI INFERRED (dashed)
  lines(seq_along(ns_inf), ns_inf, lty = 2, col = col_ns, lwd = lwd)
  
  #SSI TRUE
  lines(seq_along(ss), ss, col = col_ss, lwd = lwd)
  
  #SSI INFERRED (dashed)
  lines(seq_along(ss_inf), ss_inf, lty = 2, col = col_ss, lwd = lwd)
  
}

PLOT_SSIB_DATA_CI <- function(mcmc_ssib, list_data_ssib,
                              cex = 1.7, lwd = 2.0, 
                           col_ss = 'aquamarine', col_ns = 'orange'){
  
  #TRUE VALUES
  par(mfrow = c(1,1))
  non_ss = list_data_ssib$non_ss
  ss = list_data_ssib$ss
  
  #MCMC INFERRED
  ns_inf = mcmc_ssib$data[[1]]
  ss_inf = mcmc_ssib$data[[2]]
  lower_ci_ns <- apply(mcmc_ssib$non_ss_tot, 2, get_lower_ci)
  lower_ci_ss <- apply(mcmc_ssib$ss_tot, 2, get_lower_ci)
  upper_ci_ns <- apply(mcmc_ssib$non_ss_tot, 2, get_upper_ci)
  upper_ci_ss <- apply(mcmc_ssib$ss_tot, 2, get_upper_ci)
  mean_ci_ns <- apply(mcmc_ssib$non_ss_tot, 2, mean)
  mean_ci_ss <- apply(mcmc_ssib$ss_tot, 2, mean)
  
  #************
  #nsI TRUE
  title = bquote('Non Super-Spreading Infections - True vs Inferred')
  y_lim = c(0, max(non_ss, upper_ci_ns))
  plot(seq_along(non_ss), non_ss,  type = 'l',
       xlab = 'Time', ylab = 'Daily infection count',
       main = title, ylim = y_lim, lwd = lwd + 1.5)
       #col = col_ns, lwd = lwd)
       #cex.lab= cex, cex.axis= cex, cex.main= cex, cex.sub=cex)
  
  #nsI INFERRED (dashed)
  #Mean
  lines(seq_along(mean_ci_ns), mean_ci_ns, col = col_ns, lwd = lwd)
  #Final
  lines(seq_along(ns_inf), ns_inf, lty = 2, col = col_ns, lwd = lwd)
  
  #CIs
  x = seq_along(upper_ci_ns)
  polygon(c(x, rev(x)),
          c(upper_ci_ns, rev(lower_ci_ns)),
          col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
  
  legend('topleft', c('True', 'Mean inferred', 'Final inferred','95% CIs'),
         col = c('black', col_ns, col_ns, 'gray'), lwd = c(3, 2, 2,2), lty = c(1,1,2,1))

  #************
  #SSI TRUE
  title = bquote('Super-Spreading Infections - True vs Inferred')
  y_lim = c(0, max(ss, upper_ci_ss))
  plot(seq_along(ss), ss,  type = 'l',
       xlab = 'Time', ylab = 'Daily infection count',
       main = title, ylim = y_lim, lwd = lwd + 1.5) 
       #col = col_ss
       #cex.lab= cex, cex.axis= cex, cex.main= cex, cex.sub=cex)
  
  #SSI INFERRED (dashed)
  #Mean
  lines(seq_along(mean_ci_ss), mean_ci_ss, col = col_ss, lwd = lwd)
  #Final
  lines(seq_along(ss_inf), ss_inf, lty = 2, col = col_ss, lwd = lwd)

  #CIs
  x = seq_along(upper_ci_ss)
  polygon(c(x, rev(x)),
          c(upper_ci_ss, rev(lower_ci_ss)),
          col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
  
  legend('topleft', c('True', 'Mean', 'Final', '95% CIs'),
         col = c('black', col_ss, col_ss, 'gray'), lwd = c(3, 2, 2, 2), lty = c(1,1,2,1))
}

#PLOT
#plot.ts(mcmc_ssib$non_ss_tot[2:4800, 49], main = 'NS, Day 49', ylab = 'Daily Infection count')
