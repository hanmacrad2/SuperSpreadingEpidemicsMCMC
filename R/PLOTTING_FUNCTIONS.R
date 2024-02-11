#PLOTTING FUNCTIONS

#SIMULATED DATA 
PLOT_SIM_DATA <- function(epidemic_data, FLAGS_MODELS, RESULTS_FOLDER = '', 
                          fig_num = 1, cex = 1.9, main_font = 2.8, axis_font = 2.4,
                          PDF = FALSE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  ylabel = 'Infection count' #Daily
  
  if(PDF){
    pdf_file = paste0('/sim_', model, '_', fig_num, '.pdf') 
    pdf(paste0(RESULTS_FOLDER, model, pdf_file), width = 10.0, height = 8.0)
    
    #MARGIN
    par(mar= rep(5.0, 4), xpd=TRUE) 
  }
  
  if(FLAGS_MODELS$Baseline){
    data_title = bquote(paste(.(model), " simulated data. ", italic(R[0]),
                                   " = 2.0"))
  } else if (FLAGS_MODELS$SSE){
    
    data_title = bquote(paste(.(model), " simulated data: ", italic(R[0]),
                                   " = 2.0, ", italic(k), " = 0.1"))
  } else if (FLAGS_MODELS$SSI){
    
    data_title = bquote(paste(.(model), " simulated data: ", italic(R[0]),
                                   " = 2.0, ", k, " = 0.16"))
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
       lwd = 2.5, #3.5,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)  
  
  if(PDF){
    dev.off()
  }
}


#MCMC TRACE
PLOT_MCMC_TRACE <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                             MODEL_COLOR, cex = 1.8){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
    
  plot(seq_along(mcmc_vec), mcmc_vec,  type = 'l',
       xlab = 'Time', 
       ylab = list_labels$lab, 
       main =  list_labels$main_trace,
       col = MODEL_COLOR,
       cex.lab=cex+0.3, cex.axis=cex, cex.main=cex+0.35, cex.sub=cex)
}

#MCMC HISTOGRAM
PLOT_MCMC_HIST <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM, MODEL_COLOR,
                            cex = 1.8){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  prior = GET_PRIOR_TITLE(FLAG_PARAM)
  
  hist(mcmc_vec, freq = FALSE, breaks = 200,
       xlab = list_labels$lab,
       #xlim = xlim,
       border = MODEL_COLOR,
       col = MODEL_COLOR, 
       main = list_labels$main_hist_prior,
       cex.lab=cex+0.35, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#MCMC TRACE
PLOT_R0_TRACE <- function (r0_mcmc, FLAGS_MODELS, MODEL_COLOR, cex = 1.8){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  plot(seq_along(r0_mcmc), r0_mcmc,  type = 'l',
       xlab = 'Time', 
       ylab = expression(paste('R'[0])), 
       main =  bquote(paste(italic(R[0]), " Trace - ", .(model)~ "model")),
       col = MODEL_COLOR,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#PLOT CUMULATIVE MEAN
PLOT_CUMULATIVE_MEAN <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                             MODEL_COLOR, cex = 1.8, ylim =  c(1.8, 2.2)){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  cum_mean = cumsum(mcmc_vec)/seq_along(mcmc_vec)
  
  plot(seq_along(cum_mean), cum_mean, 
       ylim = ylim,
       col = MODEL_COLOR,
       xlab = 'Time', ylab = list_labels$lab,
       main =  list_labels$main_mean_sim,
       cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#PLOT LOG-LIKE
PLOT_LOG_LIKELIHOOD <- function(loglike_vec, FLAGS_MODELS, n_mcmc, cex = 1.8,
                                PLOT_LOGLIKE_MEAN = FALSE){
  
  'PLOT LOG LIKELIHOOD'
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  #1. TRACE
  plot(seq_along(loglike_vec), loglike_vec, type = 'l',
          xlab = 'Time', ylab = 'Log likelihood', 
          main = bquote(paste(.(model), " log likelihood. N mcmc: ", .(n_mcmc))),
          cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)
  
  if(PLOT_LOGLIKE_MEAN){
    
    #2. CUMULATIVE MEAN
    cum_mean = cumsum(loglike_vec)/seq_along(loglike_vec)
    
    plot(seq_along(cum_mean), cum_mean, 
         #ylim = ylim,
         col = 'black',
         xlab = 'Time', ylab = 'Log likelihood',
         main = bquote(paste("Log likelihood - ", .(model), " model")),
         cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
  }

}

#PLOT SIGMA
PLOT_SIGMA <- function(sigma_vec, cex = 1.8){
  
  'Plot adaptive sigma. *Need to fix in write up'
  sigma_vec = sigma_vec[2:length(sigma_vec)]
  label_sigma = bquote(paste( "Adaptive ", italic(sigma))) 
    
  plot(seq_along(sigma_vec), sigma_vec, 
       xlab = 'Time', ylab = label_sigma, 
          main =label_sigma,
          cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}

#****************************************************************************
#* MODELS

GET_FLAGS_MODELS <-function(BASELINE = FALSE, SSE = FALSE, SSI = FALSE , SSEB = FALSE, SSIB = FALSE){
  
  if(BASELINE){
    FLAGS_MODELS = list(Baseline = TRUE, SSE = FALSE, SSI = FALSE,
                        SSEB = FALSE, SSIB = FALSE) 
  } else if (SSE){
    FLAGS_MODELS = list(Baseline = FALSE, SSE = TRUE, SSI = FALSE,
                        SSEB = FALSE, SSIB = FALSE) 
  } else if (SSI) {
    
    FLAGS_MODELS = list(Baseline = FALSE, SSE = FALSE, SSI = TRUE,
                        SSEB = FALSE, SSIB = FALSE) 
    
  } else if (SSEB) {
    
    FLAGS_MODELS = list(Baseline = FALSE, SSE = FALSE, SSI = FALSE,
                        SSEB = TRUE, SSIB = FALSE) 
    
  } else if (SSIB) {
    
    FLAGS_MODELS = list(Baseline = FALSE, SSE = FALSE, SSI = FALSE,
                        SSEB = FALSE, SSIB = TRUE) 
    
  }
  
  return(FLAGS_MODELS)
}

GET_MODEL_COLORS <- function(){
  
  MODEL_COLORS <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red')
  
  return(MODEL_COLORS)
}

#PARAMS
GET_PARAM <- function(r0 = FALSE, k = FALSE, alpha = FALSE,
                            beta = FALSE, a = FALSE, b = FALSE, c = FALSE) {
  
  if(r0){
    FLAG_PARAM = list(r0 = TRUE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE)
    
  } else if (k){
    FLAG_PARAM = list(r0 = FALSE, k = TRUE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE)
    
  } else if (alpha) {
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = TRUE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE)
    
  } else if (beta){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = TRUE, a = FALSE, b = FALSE)
  } else if (a){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = TRUE, b = FALSE, c = FALSE)
  } else if (b){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = TRUE, c = FALSE)
  } else if (c) {
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = TRUE)
  }
  
  return(FLAG_PARAM)
}

#PARAM LABEL
GET_PARAM_LABEL <- function(FLAG_PARAM, model) { #model or FLAG_MODELS
  
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  if(FLAG_PARAM$r0){
    list_labels = list(lab = expression(paste('R'[0])), 
                       main_trace =  bquote(paste(italic(R[0]), " MCMC Trace")), #.(model)~ "model")),
                       main_hist = bquote(paste(italic(R[0]), " Posterior. ")),
                       main_hist_prior = bquote(paste(italic(R[0]), " Posterior. Prior: Exponential(1)")),
                       main_mean_sim = bquote(paste(italic(R[0]), " Cumulative mean. Simulated = 2.0")), #value 
                       main_mean0 = bquote(paste(italic(R[0]), " Cumulative mean - ", .(model)~ "model")))
    
  } else if (FLAG_PARAM$alpha) {
    list_labels = list(lab = expression(paste(italic(alpha))), 
                                         main_trace =  bquote(paste(italic(alpha), " MCMC Trace")), #, .(model)~ "model")),
                       main_hist = bquote(paste(italic(alpha), " Posterior. ", .(model)~ "model")),
                       main_hist_prior = bquote(paste(italic(alpha), " Posterior. Prior: Beta(2,2)")),
                       main_mean_sim = bquote(paste(italic(alpha), " Cumulative mean. Simulated = 0.5")),
                       main_mean2 = bquote(paste(italic(alpha), " Cumulative mean - ", .(model)~ "model")))
    
  } else if (FLAG_PARAM$beta){
    list_labels = list(lab = expression(paste(italic(beta))), 
                                         main_trace =  bquote(paste(italic(beta), " MCMC Trace")), #.(model)~ "model")),
                       main_hist = bquote(paste(italic(beta), " Posterior - ", .(model)~ "model")),
                       main_hist_prior = bquote(paste(italic(beta), " Posterior. Prior: 1 + Gamma(3,3)")),
                       main_mean_sim = bquote(paste(italic(beta), " Cumulative mean. Simulated = 10")),
                       main_mean2 = bquote(paste(italic(beta), " Cumulative mean - ", .(model)~ "model")))
  } else {
    
    list_labels = list(lab = bquote(paste(.(param))), 
                                         main_trace =  bquote(paste(.(param), " MCMC Trace")), #, .(model)~ "model")), #" Trace - ", .(model)~ "model")),
                       main_hist =  bquote(paste(.(param), " Posterior - ", .(model)~ "model")),
                       main_mean2 =  bquote(paste(.(param), " Cumulative mean - ", .(model)~ "model")))
    
    list_labels = GET_ADDITIONAL_TITLES(FLAG_PARAM, list_labels)
  } 
  
  return(list_labels)
}

#PRIOR TITLE FOR OTHER PARAMS
GET_ADDITIONAL_TITLES <- function(FLAG_PARAM, list_labels){
  
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  if(FLAG_PARAM$k){
    
  list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: Exponential(5)"))
  list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated = 0.1"))
  
  } else if (FLAG_PARAM$a) {
  
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: Beta(2,2)"))
    list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated = 0.5"))
    
  } else if (FLAG_PARAM$b) {
    
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: 1 + Gamma(3,3)"))
    list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated value = 10"))
  
    } else if (FLAG_PARAM$c){
      
      list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: 1 + Gamma(3,3)"))
      list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated = 10"))
  }
  
  return(list_labels)
}

#*************************
#PLOT CREDIBLE INTERVALS Nice shaded plot :) 
#*************************
# Plotting
# x = seq_along(mean_ss)
# mean_ss = unlist(as.vector(as.matrix(df_ss["mean_ss", , drop = TRUE])))
# lower_ss = unlist(as.vector(as.matrix(df_ss["lower_ci_ss", , drop = TRUE])))
# upper_ss = unlist(as.vector(as.matrix(df_ss["upper_ci_ss", , drop = TRUE])))
# 
# plot(seq_along(mean_ss), mean_ss,  pch = 16, col = "blue",
#      #ylim = range(c(df_ss$credible_lower, df_ss$credible_upper)),
#      xlab = "Time", ylab = "Values", main = "Quantile/Polygon Plot with Credible Intervals")
# 
# # Draw shaded region for credible intervals
# polygon(c(x, rev(x)),
#         c(upper_ss, rev(lower_ss)),
#         col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)