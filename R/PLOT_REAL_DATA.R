#************************
#
#PLOT REAL DATA FUNCTION
#
#************************
PLOT_MCMC_REAL_DATA <-function(epidemic_data, RESULTS_FOLDER, xlimits, data_type, main_title,
                               list_mcmc = list(Baseline = mcmc_baseline, SSE = mcmc_sse,
                                      SSI = mcmc_ssi, SSEB = mcmc_sseb, SSIB = mcmc_ssib), 
                               MODEL_COLORS, plot_margin = c(5.0, 5.2, 4.5, 1.5), cex = 2.0, 
                              PDF = TRUE,
                              plot_width = 13.5, plot_height = 12.5) {
                               #list_colors = c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', '#DC143C') ){
  
  
  'PLOT MCMC REAL DATA'
  #PDF
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, 'plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0('Real_data_', data_type,  '_mcmc_results_', '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
    
  }
  
  #SETUP
  par(mfrow = c(5, 2))
  par(mar = plot_margin)
  par(oma = c(1, 1, 5, 1)) 

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
      
    } else if (model == 'SSIB') {
      
      list_r0_vec[[model]] = mcmc_output$ssib_params_matrix[,1] 
    }
    
    else {
      
      list_r0_vec[[model]] = mcmc_output$r0_vec
    }
    
    #SSI 
    if (model == 'SSI' ) { #|| model == 'SSIB'){
      #PLOT HIST
      PLOT_HIST(list_r0_vec[[model]], model, MODEL_COLORS[i], cex, xlimits) #c(0, 4))
      
    } else if (model == 'SSIB') {
      PLOT_HIST(list_r0_vec[[model]], model, MODEL_COLORS[i], cex, c(0.25, 3.35))
    } else {
      
      #PLOT HIST
      PLOT_HIST(list_r0_vec[[model]], model, MODEL_COLORS[i], cex, xlimits)
    }

    #TRACE PLOTS
    PLOT_TRACE(list_r0_vec[[model]], model, MODEL_COLORS[i], cex)
  }
  
  #OUTER TITLE
  #main_title = bquote(main_title)
  title(main_title, outer = TRUE, cex.main = cex + 0.2) # main_font + 0.2)
  
  
  if(PDF){
    dev.off()
  }
}


#***********************************
#FUNCTIONS
#***********************************

#PLOT TRACE
PLOT_TRACE <-function(mcmc_vec, model_name, color_model, cex){
  
  #TITLE
  trace_title = bquote(paste(.(model_name) ~ "model. ", italic(R[0]),
                             " Trace"))

  ylab =  bquote(paste(R[0])) #bquote(paste(italic(R[0])))
  
  plot(seq_along(mcmc_vec), mcmc_vec,  type = 'l',
       xlab = 'Time', ylab = ylab, 
       main = trace_title,
       col = color_model,
       cex.lab= cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#PLOT SAMPLES
PLOT_HIST <-function(mcmc_vec, model_name, color_model, cex, xlimits){
  
  #TITLE
  hist_title = bquote(paste(.(model_name) ~ "model. ", italic(R[0]),
                             " Posterior Samples"))
  
  #trace_title = bquote(paste(R[0],  " Posterior; ", .(model_name) ~ "model. "))
  
  xlab =  bquote(paste(italic(R[0])))
  
  hist(mcmc_vec, freq = FALSE, breaks = 200,
       xlab = xlab,
       xlim = xlimits,
       border = color_model,
       col = color_model, 
       main = hist_title,
       cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub=cex)
  
}

#PLOT TRACE
PLOT_EPI_DATA <-function(epi_data, title, cex = 2.0){
  
  #TITLE
  trace_title = bquote(paste(.(title)))

  plot(seq_along(epi_data), epi_data,  type = 'l',
       ylab = 'Time', xlab = 'Infection count', 
       main = trace_title,
       cex.lab= cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}

#*********************************************************************************
#*
# PLOT EPIDEMIC OUTBREAK DATA + DATE
#
#****************************************************************************************************

#*******************
#PLOT_EPI_DATA_DATE_PDF
PLOT_EPI_DATA_DATE_PDF <- function(df_data, RESULTS_FOLDER, title, data_type, 
                                   plot_width = 11.5,  
                                   plot_height = 7.5, cex = 1.3, lwd_data = 2.5,
                                   X_AXIS_DATES = FALSE){
  
  #PLOT 
  plot_folder = paste0(RESULTS_FOLDER, 'plots/')
  create_folder(plot_folder)
  time_stamp = GET_CURRENT_TIME_STAMP()
  pdf_file = paste0(data_type, '_EPI_DATA_', time_stamp, '.pdf')  
  pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height)
  
  #PLOT SETTINGS
  par(oma = c(1, 1, 1, 1))
  par(mar = c(4.5,5,4,4))
  ylim = c(0, max(df_data$cases))
  
  #PLOT
  PLOT_EPI_DATA_DATE(df_data, title, cex, ylim, lwd_data = lwd_data, X_AXIS_DATES = X_AXIS_DATES)
  
  dev.off()
}

PLOT_EPI_DATA_DATE <-function(df_data, title, cex, ylim, 
                              lwd_data = 2.25, xlab = 'Report date',
                              X_AXIS_DATES = FALSE, SSE = FALSE, by_data = "3 days"){ #3 days
  
  'PLOT_EPI_DATA_DATE'
  if(X_AXIS_DATES || SSE){
    xaxt_temp = 'n'
  } else {
    xaxt_temp = 's'
  }
  
  #PLOT
  plot(df_data$date, df_data$cases, 
       type = "l", 
       xlab = xlab, 
       ylab = "Number of Daily Cases", 
       main = bquote(paste(.(title))),
       ylim = ylim, lwd = lwd_data, 
       cex.lab=cex+0.2, cex.axis=cex+0.2, cex.main=cex+0.3, cex.sub=cex+0.2,
       xaxt = xaxt_temp)
  
  if(X_AXIS_DATES){
    
    key_dates = c("2020-03-20", "2021-01-01", "2021-08-20", "2021-12-31")
    key_date_names = c("March 2020", "2021", "August 2021", "2022")
    axis.Date(side = 1, at = as.Date(key_dates), 
              labels = key_date_names, cex.axis = cex+0.2)
    
  } else if (SSE){
    dates_plot = df_data$date
    label_dates <- seq(min(dates_plot), max(dates_plot), by = by_data)
    label_dates <- label_dates[label_dates %in% dates_plot]
    axis.Date(1, at = label_dates,
             labels = format(label_dates, "%b %d"),
             cex.axis = cex + 0.2)
   }

}


#*******************
#GET_EPI_DATA_PDF
GET_EPI_DATA_PDF <- function(RESULTS_FOLDER, data_type, 
                                   plot_width = 9.5,  
                                   plot_height = 7.2){
  
  #PLOT
  plot_folder = paste0(RESULTS_FOLDER, 'plots/')
  #create_folder(plot_folder)
  time_stamp = GET_CURRENT_TIME_STAMP()
  pdf_file = paste0(data_type, '_EPI_DATA_', time_stamp, '.pdf')  
  pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height)
  
}

#*******************
#GET_PDF_SETTING
GET_PDF_SETTING <- function(RESULTS_FOLDER, data_type, result_type,
                             plot_width = 10.0,  
                             plot_height = 9.2){
  
  #PLOT 
  plot_folder = paste0(RESULTS_FOLDER, 'plots/')
  #create_folder(plot_folder)
  time_stamp = GET_CURRENT_TIME_STAMP()
  pdf_file = paste0(result_type, '_', data_type, time_stamp, '.pdf')  
  pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height)
  
}
