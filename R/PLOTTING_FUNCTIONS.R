#PLOTTING FUNCTIONS

#********************************************
#1. SIMULATED DATA 
#********************************************
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

#********************************************
#2. MCMC TRACE
#********************************************
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

#********************************************
#3. MCMC HISTOGRAM
#********************************************
PLOT_MCMC_HIST <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM, MODEL_COLOR, xlim = c(0,1),
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

#********************************************
#5. PLOT LOG-LIKELIHOOD
#********************************************
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

#********************************************
#6. PLOT SIGMA
#*********************************************
PLOT_SIGMA <- function(sigma_vec, cex = 1.8){
  
  'Plot adaptive sigma. *Need to fix in write up'
  #sigma_vec = sigma_vec[2:length(sigma_vec)]
  label_sigma = bquote(paste( "Adaptive ", italic(sigma))) 
    
  plot(seq_along(sigma_vec), sigma_vec, 
       xlab = 'Time', ylab = label_sigma, 
          main =label_sigma,
          cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.2, cex.sub=cex)
}
