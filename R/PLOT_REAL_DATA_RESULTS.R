#RESULTS REAL DATA

#PLOT TRACE
PLOT_TRACE <-function(mcmc_vec, model_name, color_model){
  
  #TITLE
  trace_title = bquote(bold(paste(.(model_name) ~ "model. ", italic(R[0]),
                                  " Trace")))
  #trace_title = bquote(bold(paste(italic(R[0]), " trace; ", .(model_name) ~ "model. ")))
  ylab =  bquote(paste(italic(R[0])))
  
  plot(seq_along(mcmc_vec), mcmc_vec,  type = 'l',
       xlab = 'Time', ylab = ylab, 
       main = trace_title,
       col = color_model,
       cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
}

#PLOT SAMPLES
PLOT_SAMPLES <-function(mcmc_vec, model_name, color_model){
  
  #TITLE
  trace_title = bquote(bold(paste(italic(R[0]),
                                  " Posterior; ", .(model_name) ~ "model. ")))
  xlab =  bquote(paste(italic(R[0])))
  hist(mcmc_vec, freq = FALSE, breaks = 200,
       xlab = xlab,
       border = color_model,
       col = color_model, 
       main = trace_title,
       cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
  
}

#PLOT CUM MEAN
PLOT_CUM_MEAN <-function(mcmc_vec, model_name, color_model){
  
  'Cumulative mean'
  
  #TITLE
  trace_title = bquote(bold(paste(italic(R[0]),
                                  " Posterior mean; ", .(model_name) ~ "")))
  ylab =  bquote(paste(italic(R[0])))
  cum_mean = cumsum(mcmc_vec)/seq_along(mcmc_vec)
  
  #ylim
  if(model_name == 'SSIB'){
    ylimits = c(0.63, 1.13)
  } else {
    ylimits = c(1.0, 1.13)
  }
  
  plot(seq_along(cum_mean), cum_mean, 
       ylim = ylimits,
       col = color_model,
       xlab = 'Time', ylab = ylab,
       main = trace_title,
       cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
}

#************************
#PLOT REAL DATA FUNCTION
#************************
PLOT_MCMC_REAL_DATA <-function(epidemic_data, OUT_FOLDER, PDF = TRUE, fig_num = '2',
                               use = 'Waitemata', list_mcmc =
                                 list(Baseline = mcmc_baseline, SSE = mcmc_sse,
                                      SSI = mcmc_ssi, SSEB = mcmc_sseb, SSIB = mcmc_ssib), 
                               list_colors = c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red') ){
  
  
  'PLOT MCMC REAL DATA'
  
  #PDF
  if(PDF){
    pdf_file = paste0('Fig_mcmc_', use, '_', fig_num, '.pdf') 
    pdf(paste0(OUT_FOLDER, 'plots/', pdf_file), width = 13.0, height = 8.0)
  }
  
  #par(mar=c(2.0, 4.0, 2.0, 2.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  #SETUP
  par(mfrow = c(3, 5))
  #Margin
  current_mar <- par("mar")
  current_mar[2] <- 5  # Change this value as needed
  par(mar = current_mar)
  
  #MODEL SPECIFIC 
  list_r0_vec = list() 
  
  #MODELS
  for (i in 1:length(list_mcmc)){
    
    model = names(list_mcmc[i]) 
    mcmc_output = list_mcmc[[model]]
    
    #MODELS - R0 VEC
    if(model == 'SSE'){
      list_r0_vec[[model]] = mcmc_output$sse_params_matrix[,1] 
      
    } else if (model == 'SSI'){
      
      list_r0_vec[[model]] = mcmc_output$ssi_params_matrix[,1] 
      
    } else {
      
      list_r0_vec[[model]] = mcmc_output$r0_vec
    }
     
    #TRACE PLOTS
    #browser()
    PLOT_TRACE(list_r0_vec[[model]], model, list_colors[i])
    
  }
  
  for (i in 1:length(list_mcmc)){
    
    model = names(list_mcmc[i]) 
    #HIST SAMPLES
    PLOT_SAMPLES(list_r0_vec[[model]], model, list_colors[i])
  }
  
  #MEAN
  for (i in 1:length(list_mcmc)){

    model = names(list_mcmc[i])
    #HIST SAMPLES
    PLOT_CUM_MEAN(list_r0_vec[[model]], model, list_colors[i])
  }
  
  if(PDF){
    dev.off()
  }
}