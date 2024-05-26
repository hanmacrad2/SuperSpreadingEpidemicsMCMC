#POSTERIOR PREDICTIVE PLOTS

#FUNCTIONS
zigzag <- function(xs) {
  sum(abs(diff(xs)))
}

#SAMPLE BASELINE
SAMPLE_BASELINE_MCMC <- function(mcmc_output, epidemic_data, n_sample_repeats = 1000,
                            PLOT = FALSE){
  
  #SAMPLE
  num_days = length(epidemic_data)
  n_mcmc = length(mcmc_output$r0_vec)
  sample_indices = sample(1:n_mcmc, n_sample_repeats)
  matrix_sim_data = matrix(nrow = n_sample_repeats, ncol = num_days)
  
  #SAMPLE PARAMETERS
  for (i in 1:n_sample_repeats){
    
    sample_index = sample_indices[i] #MCMC SAMPLE
    
    #SAMPLE PARAMETERS
    r0 = mcmc_output$r0_vec
    
    #POSTERIOR PRED DATA
    posterior_pred_data = SIMULATE_EPI_BASELINE(num_days = num_days, r0 = r0)
    
    matrix_sim_data[i, ] = posterior_pred_data
    
    #PLOT
    if(PLOT)lines(posterior_pred_data, col = 'orange')
    
  }
  
  matrix_sim_data[is.na(matrix_sim_data)] <- 0
  
  return(matrix_sim_data)
}


#SAMPLE_SSE_MCMC
SAMPLE_SSE_MCMC <- function(mcmc_output, epidemic_data, n_sample_repeats = 1000,
                            SSI = FALSE, PLOT = FALSE){
  
  #SAMPLE
  num_days = length(epidemic_data)
  
  if(SSI){
   mcmc_ss = mcmc_output$ssi_params_matrix 
  } else {
    mcmc_ss = mcmc_output$sse_params_matrix 
  }
  
  n_mcmc = length(mcmc_ss[,1])
  sample_indices = sample(1:n_mcmc, n_sample_repeats)
  matrix_sim_data = matrix(nrow = n_sample_repeats, ncol = num_days)
  
  #SAMPLE PARAMETERS
  for (i in 1:n_sample_repeats){
    
    sample_index = sample_indices[i] #MCMC SAMPLE
    
    #SAMPLE PARAMETERS 
    r0 = mcmc_ss[sample_index, 1]
    k = mcmc_ss[sample_index, 2]
    
    #POSTERIOR PRED DATA
    posterior_pred_data = SIMULATE_EPI_SSE(num_days = num_days, r0 = r0, k = k)
   
    matrix_sim_data[i, ] = posterior_pred_data
    
    #PLOT
    if(PLOT)lines(posterior_pred_data, col = 'orange')
    
  }
  
  matrix_sim_data[is.na(matrix_sim_data)] <- 0
  
  return(matrix_sim_data)
}

#SAMPLE_SSEB_MCMC
SAMPLE_SSEB_MCMC <- function(mcmc_output, epidemic_data, n_sample_repeats = 1000,
                            PLOT = FALSE){
  
  #SAMPLE
  num_days = length(epidemic_data)
  n_mcmc = length(mcmc_output$r0_vec)
  sample_indices = sample(1:n_mcmc, n_sample_repeats)
  matrix_sim_data = matrix(nrow = n_sample_repeats, ncol = num_days)
  
  
  #SAMPLE PARAMETERS
  for (i in 1:n_sample_repeats){
    
    sample_index = sample_indices[i] #MCMC SAMPLE
    
    #SAMPLE PARAMETERS
    r0 = mcmc_output$r0_vec[sample_index]
    alpha = mcmc_output$alpha[sample_index]
    beta = mcmc_output$beta[sample_index]
    
    #POSTERIOR PRED DATA
    posterior_pred_data = SIMULATE_EPI_SSEB(num_days = num_days, r0 = r0, alpha = alpha, beta = beta)
    
    matrix_sim_data[i, ] = posterior_pred_data
    
    #PLOT
    if(PLOT)lines(posterior_pred_data, col = 'orange')
    
  }
  
  matrix_sim_data[is.na(matrix_sim_data)] <- 0
  
  return(matrix_sim_data)
}

#SAMPLE_SSIB_MCMC
SAMPLE_SSIB_MCMC <- function(mcmc_output, epidemic_data, n_sample_repeats = 1000,
                             PLOT = FALSE){
  
  #SAMPLE
  num_days = length(epidemic_data)
  n_mcmc = length(mcmc_output$ssib_params_matrix[,1])
  sample_indices = sample(1:n_mcmc, n_sample_repeats)
  matrix_sim_data = matrix(nrow = n_sample_repeats, ncol = num_days)
  
  #SAMPLE PARAMETERS
  for (i in 1:n_sample_repeats){
    
    sample_index = sample_indices[i] #MCMC SAMPLE
    
    #SAMPLE PARAMETERS
    r0 = mcmc_output$ssib_params_matrix[sample_index, 1]
    a = mcmc_output$ssib_params_matrix[sample_index, 2]
    b = mcmc_output$ssib_params_matrix[sample_index, 3]
    
    #POSTERIOR PRED DATA
    posterior_pred_data = SIMULATE_EPI_SSIB(num_days = num_days, r0 = r0, a = a, b = b)
    
    matrix_sim_data[i, ] = posterior_pred_data
    
    #PLOT
    if(PLOT)lines(posterior_pred_data, col = 'orange')
    
  }
  
  matrix_sim_data[is.na(matrix_sim_data)] <- 0
  
  return(matrix_sim_data)
}
#POSTERIOR_PREDICTIVE_PLOTS
POSTERIOR_PREDICTIVE_PLOTS <- function(matrix_sim_data, true_data, df_data,
                                       title, data_type, FLAGS_MODELS, 
                                         MODEL_COLOR, RESULTS_FOLDER, EPI_DATA = TRUE,
                                       PDF = TRUE, ZIG_ZAG = FALSE, ALL_MODELS = FALSE,
                                         plot_width = 9.0, plot_height = 6.5,
                                       cex = 1.0){
  
    #SETUP
    num_days = length(true_data)
    model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
    n_samples = nrow(matrix_sim_data)
    
    #TEST STATS
    mean_est <- apply(matrix_sim_data, 2, mean) |> unlist()
    upper_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.975) |> unlist()
    lower_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.025) |> unlist()
    posterior_zig_zag = apply(matrix_sim_data, 1, zigzag)
    zigzag_true = zigzag(true_data)
    
    #PDF
    if(PDF){
      #PLOT
      print('PDF TRUE')
      plot_folder = paste0(RESULTS_FOLDER, 'plots/')
      create_folder(plot_folder)
      time_stamp = GET_CURRENT_TIME_STAMP()
      pdf_file = paste0('POSTERIOR_PRED_PLOT_', model, '_', time_stamp, '.pdf')  
      pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height)
      
      #PLOT SETTINGS
      par(oma = c(1, 1, 1, 1))
      par(mar = c(4.5,5,4,4))
    }
    
    #PLOTS 
    print('upper_bounds'); print(upper_bounds)
    print('lower_bounds'); print(lower_bounds)
   
    #PLOT 95 % Credible intervals
    if(ZIG_ZAG){
      par(mfrow = c(1, 2)) 
    }
    
    #PLOT
    main_title = bquote(paste('Posterior Predictive Plots:', .(title)))
    ylim = c(0, max(true_data, upper_bounds))
    
    if(EPI_DATA){
      plot(1:num_days, true_data, type = 'l', ylim = ylim,
           main = main_title, xlab = 'Time', ylab = 'Infection count', lwd = 2,
           cex.lab=cex+0.2, cex.axis=cex+0.2, cex.main=cex+0.3, cex.sub=cex+0.2)
    } else {
    }
         
    #PLOT EPI DATA DATE
    PLOT_EPI_DATA_DATE(df_data, main_title, cex, ylim = ylim)
    
    #95 % Credible intervals 
    if (EPI_DATA){
      lines(1:num_days, upper_bounds, col = MODEL_COLOR, lwd = 2.5) #, type = 'l') #, ylim = c(-5, 5))
      lines(1:num_days, lower_bounds, col = MODEL_COLOR, lwd = 2.5) #, type = 'l')
    } else {
      lines(df_data$date, upper_bounds, col = MODEL_COLOR, lwd = 2.5) #, type = 'l') #, ylim = c(-5, 5))
      lines(df_data$date, lower_bounds, col = MODEL_COLOR, lwd = 2.5) #, type = 'l')
    }
    
    #LEGEND
    if(ALL_MODELS){
      GET_LEGEND_POST_PRED_ALL_MODELS(data_type, n_samples)
    } else {
      GET_LEGEND_POST_PRED(FLAGS_MODELS, MODEL_COLOR, data_type, n_samples)
    }
    
    
    #HISTOGRAM OF ZIG-ZAG 
    if(ZIG_ZAG){
      max_x_lim = max(max(zigzag_true), max(posterior_zig_zag)); min_x_lim =  min(min(zigzag_true), min(posterior_zig_zag))
      hist(posterior_zig_zag, breaks = 50, xlim = c(min_x_lim, max_x_lim), #xlim = c(0, max(max(posterior_zig_zag), zigzag(true_data) * 2)),
           main = 'Zig-zag of data: sum(abs(diff(data)))')
      abline(v = quantile(posterior_zig_zag, probs = c(0.025, 0.975)), col = 'red', lwd = 2)
      abline(v = zigzag(true_data), col = 'green', lwd = 2) 
    }
    
    if(PDF && ALL_MODELS == FALSE){
      print('PDF && ALL_MODELS == FALSE')
      dev.off()
    }
  }
  
  
#LEGEND  
GET_LEGEND_POST_PRED <- function(FLAGS_MODELS, MODEL_COLOR, data_type, n_samples,
                                 legend_location = 'topleft', 
                                 cex = 1.0, inset = 0.07){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  #PARAM
  legend_caption = c(data_type, paste0(' 95 % CIs of ', model, ' model fits, N = ', (n_samples)))
  #legend_caption = GET_OFFSPRING_LEGEND_CAPTION(list_params, FLAGS_MODELS)

  #Legend
  legend(legend_location,
         x.intersp = 0.05, #-0.05,#CONTROLS RELATIVE TO THE PLOT
         legend_caption,
         cex = cex,
         #inset = c(inset,0),
         #inset = c(-inset), #CONTROLS RELATIVE TO THE MARGINS NOT PLOT!
         col = c('black', MODEL_COLOR),
         lwd = 1.5, #rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = 1, # rep(1, num_conds), #c(1, 1),
         text.font = 1.5, #1.45
         bty = "n")
}

#LEGEND  
GET_LEGEND_POST_PRED_ALL_MODELS <- function(data_type, n_samples,
                                 legend_location = 'topleft', 
                                 models = c('Baseline', 'SSE', 'SSI', 'SSEB', 'SSIB'),
                                 MODEL_COLORS = c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', '#DC143C'), 
                                 cex = 1.15, inset = 0.07){
  
  #PARAM
  legend_caption = c(data_type)
  
  for (i in 1:length(MODEL_COLORS)){
    model = models[i]
    #model = names(FLAGS_MODELS)[unlist(FLAGS_MODELS)[[1]]][i] 
    print(model)
    legend_caption = c(legend_caption, paste0(' 95 % CIs of ', model, ' model fits, N = ', (n_samples)))
  }
  
  #Legend
  legend(legend_location,
         x.intersp = 0.05, #-0.05,#CONTROLS RELATIVE TO THE PLOT
         legend_caption,
         cex = cex,
         #inset = c(inset,0),
         #inset = c(-inset), #CONTROLS RELATIVE TO THE MARGINS NOT PLOT!
         col = c('black', MODEL_COLORS),
         lwd = 2.0, #rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = 1.5, # rep(1, num_conds), #c(1, 1),
         text.font = 1.2, #1.45
         bty = "n")
}

#GET CI DATA
GET_CI_PRED_DATA <- function(matrix_sim_data, df_data, MODEL_COLOR){
  
  #TEST STATS
  mean_est <- apply(matrix_sim_data, 2, mean) |> unlist()
  upper_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.975) |> unlist()
  lower_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.025) |> unlist()
  
  #PLOT
  lines(df_data$date, upper_bounds, col = MODEL_COLOR, lwd = 2.5) 
  lines(df_data$date, lower_bounds, col = MODEL_COLOR, lwd = 2.5) 
}


