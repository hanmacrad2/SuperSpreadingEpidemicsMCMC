#PLOTTING FUNCTIONS

#SIMULATED DATA 
PLOT_SIM_DATA <- function(epidemic_data, FLAGS_MODELS, RESULTS_FOLDER = '', 
                          fig_num = 1, cex = 1.6, main_font = 2.8, axis_font = 2.4,
                          PDF = FALSE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  ylabel = 'Daily Infection count'
  
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
    
    data_title = bquote(paste(.(model), " simulated data. ", italic(R[0]),
                                   " = 2.0, ", italic(k), " = 0.16"))
  } else if (FLAGS_MODELS$SSI){
    
    data_title = bquote(bold(paste(.(model), " sim data. ", italic(R[0]),
                                   " = 2.0, ", k, " = 0.15")))
  } else if (FLAGS_MODELS$SSEB){
    data_title = bquote(bold(paste(.(model), " sim data. ", italic(R[0]),
                                   " = 2.0, ", alpha, " = 0.5, ", beta, " = 10")))
    main_font = 2.7
    
  } else if (FLAGS_MODELS$SSIB){
    data_title = bquote(bold(paste(.(model), "sim data. ", italic(R[0]),
                                   " = 2.0, ", a, " = 0.5, ", b, " = 10")))
    main_font = 2.7
  }
  
  #PLOT
  plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
       xlab = 'Time', ylab = ylabel,
       main = data_title,
       lwd = 2.5, #3.5,
       cex.lab= cex, cex.axis= cex, cex.main= cex, cex.sub=cex)
  
  if(PDF){
    dev.off()
  }
}


#MCMC TRACE
PLOT_MCMC_TRACE <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                             MODEL_COLOR, cex = 1.6){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
    
  plot(seq_along(mcmc_vec), mcmc_vec,  type = 'l',
       xlab = 'Time', 
       ylab = list_labels$lab, 
       main =  list_labels$main_trace,
       col = MODEL_COLOR,
       cex.lab= cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#MCMC HISTOGRAM
PLOT_MCMC_HIST <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM, MODEL_COLOR,
                            cex = 1.6, xlim = c(1.7, 2.3)){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  
  hist(mcmc_vec, freq = FALSE, breaks = 200,
       xlab = list_labels$lab,
       xlim = xlim,
       border = MODEL_COLOR,
       col = MODEL_COLOR, 
       main = list_labels$main_hist,
       cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#MCMC TRACE
PLOT_R0_TRACE <- function (r0_mcmc, FLAGS_MODELS, MODEL_COLOR, cex = 1.6){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  plot(seq_along(r0_mcmc), r0_mcmc,  type = 'l',
       xlab = 'Time', 
       ylab = expression(paste('R'[0])), 
       main =  bquote(paste(italic(R[0]), " Trace - ", .(model)~ "model")),
       col = MODEL_COLOR,
       cex.lab= cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

PLOT_HISTOGRAM <-function(mcmc_vec, FLAGS_MODELS, MODEL_COLOR, cex = 1.6){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  #TITLE
  title = bquote(paste(italic(R[0]),
                                  " Posterior - ", .(model) ~ "model"))
  
  xlab =  bquote(paste(italic(R[0])))
  hist(mcmc_vec, freq = FALSE, breaks = 200,
       xlab = xlab,
       xlim = c(1.8, 2.2),
       border = MODEL_COLOR,
       col = MODEL_COLOR, 
       main = title,
       cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
  
}

#PLOT CUMULATIVE MEAN
PLOT_CUMULATIVE_MEAN <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM,
                             MODEL_COLOR, cex = 1.6, ylim =  c(1.8, 2.2)){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  cum_mean = cumsum(mcmc_vec)/seq_along(mcmc_vec)
  
  plot(seq_along(cum_mean), cum_mean, 
       ylim = ylim,
       col = MODEL_COLOR,
       xlab = 'Time', ylab = list_labels$lab,
       main =  list_labels$main_mean,
       cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

PLOT_CUMULATIVE_MEAN0 <-function(mcmc_vec, FLAGS_MODELS, MODEL_COLOR, cex = 1.6){
  
  'Cumulative mean'
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  ylab =  bquote(paste(italic(R[0])))
  cum_mean = cumsum(mcmc_vec)/seq_along(mcmc_vec)
  ylimits = c(1.8, 2.2) 
  
  plot(seq_along(cum_mean), cum_mean, 
       ylim = ylimits,
       col = MODEL_COLOR,
       xlab = 'Time', ylab = ylab,
       main =  bquote(paste(italic(R[0]) ~ " Cumulative mean ")),
       cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

PLOT_CUMULATIVE_MEAN_K <-function(mcmc_vec, FLAGS_MODELS, MODEL_COLOR, cex = 1.6){
  
  'Cumulative mean'
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  ylab =  bquote(paste(italic(k)))
  cum_mean = cumsum(mcmc_vec)/seq_along(mcmc_vec)
  #ylimits = c(1.8, 2.2) 
  
  plot(seq_along(cum_mean), cum_mean, 
       #ylim = ylimits,
       col = MODEL_COLOR,
       xlab = 'Time', ylab = ylab,
       main =  bquote(paste(italic(k) ~ " Cumulative mean ")),
       cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#PLOT LOG-LIKE
PLOT_LOG_LIKELIHOOD <- function(loglike_vec, FLAGS_MODELS, cex = 1.6){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  plot(seq_along(loglike_vec), loglike_vec, type = 'l',
          xlab = 'Time', ylab = 'Log likelihood', 
          main = bquote(paste("Log likelihood - ", .(model), " model")),
          cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#PLOT SIGMA
PLOT_SIGMA <- function(sigma_vec, cex = 1.6){
  
  sigma_vec = sigma_vec[2:length(sigma_vec)]
  
  plot(seq_along(sigma_vec), sigma_vec, 
       xlab = 'Time', ylab = 'sigma', 
          main = bquote(paste("Adaptive sigma")),
          cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#*******************
#* MODELS

GET_FLAGS_MODELS <-function(BASELINE = FALSE, SSE = FALSE, SSI = FALSE , SSEB = FALSE, SSIB = FALSE){
  
  if(BASELINE){
    FLAGS_MODELS = list(Baseline = TRUE, SSE = FALSE, SSI = FALSE,
                        SSEB = FALSE, SSIB = FALSE) 
  } else if (SSE){
    FLAGS_MODELS = list(Baseline = FALSE, SSE = TRUE, SSI = FALSE,
                        SSEB = FALSE, SSIB = FALSE) 
  }
  
  return(FLAGS_MODELS)
}

GET_MODEL_COLORS <- function(){
  
  MODEL_COLORS <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red')
  
  return(MODEL_COLORS)
}

#PARAMS
GET_PARAM <- function(r0 = FALSE, k = FALSE, alpha = FALSE,
                            beta = FALSE, a = FALSE, b = FALSE) {
  
  if(r0){
    FLAG_PARAM = list(r0 = TRUE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE)
    
  } else if (k){
    FLAG_PARAM = list(r0 = FALSE, k = TRUE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE)
    
  } else if (alpha) {
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = TRUE,
                      beta = FALSE, a = FALSE, b = FALSE)
    
  } else if (beta){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = TRUE, a = FALSE, b = FALSE)
  } else if (a){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = TRUE, b = FALSE)
  } else if (b){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = TRUE)
  }
  
  return(FLAG_PARAM)
}

#PARAM LABEL
GET_PARAM_LABEL <- function(FLAG_PARAM, model) { #model or FLAG_MODELS
  
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  if(FLAG_PARAM$r0){
    list_labels = list(lab = expression(paste('R'[0])), 
                       main_trace =  bquote(paste(italic(R[0]), " MCMC Trace - ", .(model)~ "model")),
                       main_hist = bquote(paste(italic(R[0]), " Posterior - ", .(model)~ "model")),
                       main_mean = bquote(paste(italic(R[0]), " Cumulative mean - ", .(model)~ "model")))
    
  } else if (FLAG_PARAM$alpha) {
    list_labels = list(lab = expression(paste(italic(alpha))), 
                                         main_trace =  bquote(paste(italic(alpha), " MCMC Trace - ", .(model)~ "model")),
                       main_hist = bquote(paste(italic(alpha), " Posterior - ", .(model)~ "model")),
                       main_mean = bquote(paste(italic(alpha), " Cumulative mean - ", .(model)~ "model")))
    
  } else if (FLAG_PARAM$beta){
    list_labels = list(lab = expression(paste(italic(beta))), 
                                         main_trace =  bquote(paste(italic(beta), " MCMC Trace - ", .(model)~ "model")),
                       main_hist = bquote(paste(italic(beta), " Posterior - ", .(model)~ "model")),
                       main_mean = bquote(paste(italic(beta), " Cumulative mean - ", .(model)~ "model")))
  } else {
    
    list_labels = list(lab = bquote(paste(.(param))), 
                                         main_trace =  bquote(paste(.(param), " Trace - ", .(model)~ "model")),
                       main_hist =  bquote(paste(.(param), " Posterior - ", .(model)~ "model")),
                       main_mean =  bquote(paste(.(param), " Cumulative mean - ", .(model)~ "model")))
  } 
  
  return(list_labels)
}


