#PLOTTING FUNCTIONS

#********************************************
#1. SIMULATED DATA 
#********************************************
PLOT_SIM_DATA <- function(epidemic_data, FLAGS_MODELS, RESULTS_FOLDER,
                          cex = 1.1, main_font = 2.5, axis_font = 1.6,
                          PDF = TRUE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  ylabel = 'Infection count' #Daily
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_EPI_DATA_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = 7.5, height = 5.2) 
  }
  
  par(mar = c(5, 5, 4, 1.5)) #bottom, left, top, right
  
  if(FLAGS_MODELS$Baseline){
    data_title = bquote(paste(.(model), " simulated epidemic data. ", italic(R[0]),
                              " = 2.0"))
  } else if (FLAGS_MODELS$SSE){
    
    data_title = bquote(paste(.(model), " simulated epidemic data: ", italic(R[0]),
                              " = 2.0, ", italic(k), " = 0.1"))
  } else if (FLAGS_MODELS$SSI){
    
    data_title = bquote(paste(.(model), " simulated epidemic data: ", italic(R[0]),
                              " = 2.0, ",italic(k), " = 0.2 "))
  } else if (FLAGS_MODELS$SSEB){
    data_title = bquote(paste(.(model), " simulated data: ", italic(R[0]),
                              " = 2.0, ", alpha, " = 0.5, ", beta, " = 10"))
    main_font = 2.7
    
  } else if (FLAGS_MODELS$SSIB){
    data_title = bquote(paste(.(model), " simulated data: ", italic(R[0]),
                              " = 2.0, ", a, " = 0.5, ", b, " = 10"))
    main_font = 2.7
  }
  
  #PLOT
  plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
       xlab = 'Time', ylab = ylabel,
       main = data_title,
       lwd = 2.0, #3.5,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)  
  
  if(PDF){
    dev.off()
  }
}

#********************************************
#2. MCMC TRACE
#********************************************
PLOT_MCMC_TRACE <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                             RESULTS_FOLDER, MODEL_COLOR, sim_val,
                             cex = 1.1, PDF = TRUE, FLAG_DATA_AUG = FALSE){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  
  if(FLAG_DATA_AUG){
    FLAG_MODEL_DA = GET_FLAGS_MODELS_2(SSIB_NO_DA = TRUE)
    list_labels = GET_PARAM_LABEL_DA(FLAG_PARAM, FLAG_MODEL, FLAG_MODEL_DA) 
  }
  
  #PLOTTING
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_', param, '_MCMC_TRACE_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = 7.5, height = 5.2) 
  }
  par(mar = c(5, 5, 4, 1.5)) #bottom, left, top, right
  
  plot(seq_along(mcmc_vec), mcmc_vec,  type = 'l',
       xlab = 'Time', 
       ylab = list_labels$lab, 
       main =  list_labels$title_trace, #list_labels$main_trace,
       col = MODEL_COLOR,
       cex.lab=cex+0.3, cex.axis=cex, cex.main=cex+0.35, cex.sub=cex)
  
  #PLOT_TRUE_LINE
  PLOT_TRUE_LINE(sim_val, mcmc_vec) 
  
  if(PDF){
   dev.off() 
  }
}

#********************************************
#3. MCMC HISTOGRAM
#********************************************
PLOT_MCMC_HIST <- function (mcmc_vec, list_labels, FLAG_PARAM,
                            MODEL_COLOR, xlim, ylim, 
                            cex = 1.4, breaks_hist = 100, single_inf = TRUE){
  
  prior = GET_PRIOR_TITLE(FLAG_PARAM)
  
  if(single_inf){
    COLOR_BORDER = 'black'
    hist_title = list_labels$title_hist_inf
  } else {
    COLOR_BORDER = MODEL_COLOR
    hist_title = list_labels$main_hist_prior
  }
  
  hist(mcmc_vec, freq = FALSE, breaks = breaks_hist,
       xlab = list_labels$lab,
       xlim = xlim,
       ylim = ylim,
       border = COLOR_BORDER,
       col = MODEL_COLOR, 
       main = hist_title,
       cex.lab=cex+0.35, cex.axis=cex+0.3, cex.main=cex+0.2, cex.sub=cex+0.3)
}

#*************************
#* HIST + PRIOR 
PLOT_MCMC_HIST_AND_PRIOR <-function(mcmc_vec, FLAG_MODEL, FLAG_PARAM, RESULTS_FOLDER,
                                    MODEL_COLOR, sim_val, xlim, ylim,
                                    inset =  -0.007, cex = 1.75,
                                    breaks_hist = 100, legend_location = 'topright', 
                                    diff_max_y_bar_y = 0.05, FLAG_DATA_AUG = FALSE,
                                    PDF = TRUE ){
  #PARAMS
  model = names(FLAG_MODEL)[which(unlist(FLAG_MODEL))] 
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  prior_title = GET_PRIOR_TITLE(FLAG_PARAM)
  legend_list = c(list_labels$legend_mcmc_hist, prior_title,
                  list_labels$legend_true_val)
  true_bar_height = ylim[2] - diff_max_y_bar_y
  
  if(FLAG_DATA_AUG){
    FLAG_MODEL_DA = GET_FLAGS_MODELS_2(SSIB_NO_DA = TRUE)
    list_labels = GET_PARAM_LABEL_DA(FLAG_PARAM, FLAG_MODEL, FLAG_MODEL_DA) 
  }
  
  #PLOTTING
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_', param, '_HISTOGRAM_', time_stamp, '.pdf')  
    #pdf(paste0(plot_folder, pdf_file), width = 11.0, height = 7.2) 
    #pdf(paste0(plot_folder, pdf_file), width = 11.5, height = 7.5) 
    pdf(paste0(plot_folder, pdf_file), width = 14.5, height = 9.0) 
  }
  #MARGIN
  par(mar = c(5, 5, 4, 2)) #bottom, left, top, right #1.5
  
  #HIST
  PLOT_MCMC_HIST(mcmc_vec, list_labels, FLAG_PARAM, 
                 MODEL_COLOR, xlim, ylim, cex = cex, breaks_hist = breaks_hist)
  
  #PRIOR
  PLOT_PRIOR_DIST(FLAG_PARAM, xlim, alpha = 0.4)
  
  #PLOT TRUE VERTICAL LINE 
  PLOT_TRUE_LINE(sim_val, mcmc_vec, true_bar_height, VERTICAL = TRUE)
  
  #LEGEND
  GET_LEGEND_INF_HIST(legend_list, MODEL_COLOR, legend_location, inset)
  
  if(PDF){
    dev.off()
  }
  
}

#********************************************
#3b. R0 MCMC TRACE
#********************************************
PLOT_R0_TRACE <- function (r0_mcmc, FLAGS_MODELS, MODEL_COLOR, cex = 1.8){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  plot(seq_along(r0_mcmc), r0_mcmc,  type = 'l',
       xlab = 'Time', 
       ylab = expression(paste('R'[0])), 
       main =  bquote(paste(italic(R[0]), " Trace - ", .(model)~ "model")),
       col = MODEL_COLOR,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#********************************************
#4. PLOT CUMULATIVE MEAN
#********************************************
PLOT_CUMULATIVE_MEAN <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                                  MODEL_COLOR, cex = 1.8, ylim =  c(1.8, 2.2), SINGLE_INF = TRUE){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  cum_mean = cumsum(mcmc_vec)/seq_along(mcmc_vec)
  
  if(SINGLE_INF){
    main_title = list_labels$cum_mean
  } else {
    main_title = list_labels$main_mean_sim
  }
  
  plot(seq_along(cum_mean), cum_mean, 
       ylim = ylim,
       col = MODEL_COLOR,
       xlab = 'Time', ylab = list_labels$lab,
       main =  main_title,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#********************************************
#5. PLOT LOG-LIKELIHOOD
#********************************************
PLOT_LOG_LIKELIHOOD <- function(loglike_vec, FLAGS_MODELS, RESULTS_FOLDER,
                                COL_LOG_LIKE,
                                cex = 1.1, n_mcmc = 125000,
                                #COL_LOG_LIKE = 'black',
                                PDF = TRUE, PLOT_LOGLIKE_MEAN = FALSE,
                                FLAG_DATA_AUG = FALSE){
  
  'PLOT LOG LIKELIHOOD'
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  #TITLE OF PLOT
  if(FLAG_DATA_AUG){
    title = bquote(paste("The log-likelihood of the ", .(model), " model, Known SS Infections"))
  } else {
    title = bquote(paste("The log-likelihood of the ", .(model), " model"))
  }
  
  #PLOTTING
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_LOG_LIKE_TRACE_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = 7.5, height = 5.2) 
  }
  par(mar = c(5, 5, 4, 1.5)) #bottom, left, top, right
  
  #1. TRACE
  plot(seq_along(loglike_vec), loglike_vec, type = 'l',
       xlab = 'Time', ylab = 'Log likelihood', 
       col = COL_LOG_LIKE,
       main = title,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)
  
  if(PDF){
    dev.off()
  }
  if(PLOT_LOGLIKE_MEAN){
    
    #2. CUMULATIVE MEAN
    cum_mean = cumsum(loglike_vec)/seq_along(loglike_vec)
    
    plot(seq_along(cum_mean), cum_mean, 
         #ylim = ylim,
         col = 'black',
         xlab = 'Time', ylab = 'Log likelihood',
         main = bquote(paste("Log likelihood - ", .(model), " model")), 
         #main = bquote(paste(.(model), " log likelihood. N mcmc: ", .(n_mcmc))),
         cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
  }
  
}

#********************************************
#6. PLOT SIGMA
#*********************************************
PLOT_SIGMA <- function(sigma_vec, cex = 1.8){
  
  'Plot adaptive sigma. *Need to fix in write up'
  sigma_vec = sigma_vec[2:length(sigma_vec)]
  label_sigma = bquote(paste( "Adaptive ", italic(sigma))) 
  
  plot(seq_along(sigma_vec), sigma_vec, 
       xlab = 'Time', ylab = label_sigma, 
       main =label_sigma,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#************
#GET LEGEND
#************
GET_LEGEND_INF_HIST <- function(legend_list, MODEL_COLOR, # inset = 0.25 # legend_location = 'topright',  inset = 0.5)
                                   legend_location, inset){
  
  #COLOUR
  COLOR_ALPHA = GET_COLOR_ALPHA(MODEL_COLOR, alpha = 0.8)
  
  #Legend
  num_conds = length(legend_list)
  pch_list = rep(19, num_conds)
  legend(legend_location, #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         cex = 1.45,
         inset=c(inset,0), #c(-inset,0),
         col = c(COLOR_ALPHA, 'gray', 'black'),
         lwd = rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = rep(1, num_conds), #c(1, 1),
         #pch = pch_list, #c(NA, pch_list, NA, NA, NA),
         text.font = 2.3, #1.45
         pt.cex = 0.7,
         bty = "n")
}

#********************************************
#3. MCMC HISTOGRAM
#********************************************
PLOT_MCMC_HIST_V0 <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                               MODEL_COLOR, xlim, ylim, 
                               cex = 1.4, breaks_hist = 100, single_inf = TRUE){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  prior = GET_PRIOR_TITLE(FLAG_PARAM)
  
  if(single_inf){
    COLOR_BORDER = 'black'
    hist_title = list_labels$title_hist_inf
  } else {
    COLOR_BORDER = MODEL_COLOR
    hist_title = list_labels$main_hist_prior
  }
  
  hist(mcmc_vec, freq = FALSE, breaks = breaks_hist,
       xlab = list_labels$lab,
       xlim = xlim,
       ylim = ylim,
       border = COLOR_BORDER,
       col = MODEL_COLOR, 
       main = hist_title,
       cex.lab=cex+0.35, cex.axis=cex+0.3, cex.main=cex+0.2, cex.sub=cex+0.3)
}

#SSIB MODEL
PLOT_SSIB_DATA_CI <- function(mcmc_ssib, list_data_ssib, RESULTS_FOLDER,
                              cex = 1.2, lwd = 2.0, model = 'SSIB',
                              col_ss = 'aquamarine', col_ns = 'orange', PDF = TRUE){
  
  #PLOT
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0('NS_data_', model, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = 7.75, height = 5.6) 
  }
  par(mar = c(5, 5, 4, 2)) #bottom, left, top, right #1.5
  
  #TRUE VALUES
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
  mean_ci_ns <- apply(mcmc_ssib$non_ss_tot, 2, median)
  mean_ci_ss <- apply(mcmc_ssib$ss_tot, 2, median)
  
  #************
  #NS TRUE
  title = bquote(paste('Non Super-Spreading Infections, ', .(model), ' model - True vs Inferred'))
  y_lim = c(0, max(non_ss, upper_ci_ns))
  plot(seq_along(non_ss), non_ss,  type = 'l',
       xlab = 'Time', ylab = 'Infection count',
       main = title, 
       ylim = y_lim, 
       lwd = 2.0,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)  
  
  #NS INFERRED (dashed)
  #Mean
  lines(seq_along(mean_ci_ns), mean_ci_ns, col = col_ns, lwd = lwd)
  
  #Final DOTTED
  lines(seq_along(ns_inf), ns_inf, lty = 2, col = col_ns, lwd = lwd)
  
  #CIs
  x = seq_along(upper_ci_ns)
  polygon(c(x, rev(x)),
          c(upper_ci_ns, rev(lower_ci_ns)),
          col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
  
  #mean
  legend('topleft', c('True value', 'Median inferred',
                      'Final inferred','95% CIs'),
         col = c('black', col_ns, col_ns, 'gray'), lwd = c(3, 2, 2,2), lty = c(1,1,2,1))
  
  #************
  #SSI TRUE
  
  if(PDF){
    dev.off()
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0('SS_data_', model, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), 7.75, height = 5.6) 
  }
  par(mar = c(5, 5, 4, 2)) #bottom, left, top, right #1.5
  
  title = bquote(paste('Super-Spreading Infections, ', .(model), ' model - True vs Inferred'))
  #title = bquote('Super-Spreading Infections - True vs Inferred')
  y_lim = c(0, max(ss, upper_ci_ss))
  plot(seq_along(ss), ss,  type = 'l',
       xlab = 'Time', ylab = 'Daily infection count',
       main = title, 
       ylim = y_lim,
       lwd = 2.0,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)  
  
  #SSI INFERRED (dashed)
  #Mean
  lines(seq_along(mean_ci_ss), mean_ci_ss, col = col_ss, lwd = lwd)
  
  #Final (DOTTED)
  lines(seq_along(ss_inf), ss_inf, lty = 2, col = col_ss, lwd = lwd)
  
  #CIs
  x = seq_along(upper_ci_ss)
  polygon(c(x, rev(x)),
          c(upper_ci_ss, rev(lower_ci_ss)),
          col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
  
  legend('topleft', c('True value', 'Median inferred', 'Final inferred', '95% CIs'),
         col = c('black', col_ss, col_ss, 'gray'), lwd = c(3, 2, 2, 2),
         lty = c(1,1,2,1))
  
  if(PDF){
    dev.off()
  }
}

#PLOT LINE TRUE VALUE
PLOT_TRUE_LINE <- function(sim_val, mcmc_vec, true_bar_height, VERTICAL = FALSE){
  
  if(VERTICAL){
  segments(sim_val, 0, sim_val, true_bar_height, col = 'black', lwd = 3.0)
  } else {
    segments(0, sim_val, length(mcmc_vec), sim_val, col = 'black', lwd = 2)
  }
} 

