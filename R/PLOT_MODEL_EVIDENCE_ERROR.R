#*********************************
#*
#* PLOT MODEL EVIDENCE ERROR
#* 
#* ******************************

PLOT_VIOLIN_MODEL_EVIDENCE <- function(vec_results, FLAGS_MODELS, 
                                         MODEL_COLOR, RESULTS_FOLDER,
                                       FLAG_RESULT_TYPE = list('Importance Sampling' = FALSE, 
                                                               'Harmonic Mean' = TRUE),
                                       y_lim = c(-131, -129),
                                         plot_width = 6.5, plot_height = 5.5,
                                         cex = 1.9, main_font = 2.2, #3.2, 
                                         axis_font = 1.6, PDF = TRUE){
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  result_type = names(FLAG_RESULT_TYPE)[which(unlist(FLAG_RESULT_TYPE))]
  main_title =  bquote(paste(.(model), " model - ", .(result_type), " Estimated Evidence"))
  #param = list_params[1]
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_', result_type, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  #PLOT
  if(FLAGS_MODELS$SSIB){
    par(mar = c(4.5, 4.5, 3, 2)) #c(5, 4.7, 4, 4.7))
  } else {
    par(mar = rep(5,4))
  }
  
  #par(mar = c(5.5, 4.7, 4, 4.7))
  #par(oma = c(1, 1, 1, 1)) #bottom, left, top, right
  
  violin_plot(vec_results, #names.arg = x, 
          col = MODEL_COLOR, #"skyblue",
          main = main_title,
          #xlab = bquote(paste(.(model), " model")),
          ylab = "log Model Evidence",
          #lwd = lwd, #3.5,
          ylim = y_lim,
          cex.lab=cex+0.35, cex.axis=cex + 0.2, 
          cex.main=cex+0.3, cex.main= main_font,
          x_axis_labels = c(""), )
         # xaxt = "n") 
  
  GET_LEGEND_SD_ERROR(sd(vec_results), MODEL_COLOR)
  
  mtext(side = 1, line = 2, 
        bquote(paste(.(model), " model")),
        cex = 1.3)
  
  #axis(side = 1, at = 1, labels = "", tick = FALSE)
  
  #GET_MODEL_PLOT_SETTINGS(FLAGS_MODELS, FLAGS_PARAM)
  
  #MAIN TITLE
  #title(main_title, outer = TRUE, cex.main = main_font + 0.2)
  
  if(PDF){
    dev.off()
  }
}

#GET LEGEND SD OR THE ERROR
GET_LEGEND_SD_ERROR <- function(sd_error, MODEL_COLOR,
                                 legend_location = 'bottomright', #,'topright', 
                                 cex = 1.1, inset = 0.1){
  
  #Legend
  legend_caption = c('N repeats = 1000',
                     paste0('sd = ', round(sd_error, 3)))
  
  legend(legend_location,
         #x.intersp = 0.05, #-0.05,#CONTROLS RELATIVE TO THE PLOT
         legend_caption,
         cex = cex,
         #inset = -0.02,
         #inset = c(-inset), #CONTROLS RELATIVE TO THE MARGINS NOT PLOT!
         col = c(MODEL_COLOR, 'white'),
         lwd = c(2,2), #rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = c(1,1), # rep(1, num_conds), #c(1, 1),
         text.font = 1.5, #1.45
         bty = "n")
}
