#********************************************
#1. PLOT_OFF SPRING DISTRIBUTIONS
#********************************************
PLOT_SIMULATION_LIST <- function(list_epidemic_data, list_params,
                                 FLAGS_MODELS, FLAGS_PARAM, 
                                 MODEL_COLOR, RESULTS_FOLDER,
                                 #plot_width = 11.0, plot_height = 9.0,
                                 plot_width = 13.5, plot_height = 11.0,
                                 cex = 2.3, lwd = 2.45,
                                 main_font = 3.5, #3.2, 
                                 axis_font = 1.6,
                                 PDF = TRUE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAGS_PARAM)[which(unlist(FLAGS_PARAM))]
  ylabel = 'Infection count' #Daily
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_EPI_DATA_', param, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  #PLOT SETTINGS
  GET_MODEL_PLOT_SETTINGS(FLAGS_MODELS, FLAGS_PARAM)
  
  #PLOT
  for(i in 1:length(list_epidemic_data)){
    
    #PARAMS
    list_titles = GET_MODEL_SIMULATION_TITLE(i, list_params, 
                                             FLAGS_MODELS, FLAGS_PARAM)
    data_title = list_titles$data_title
    main_title = list_titles$main_title
    
    epidemic_data = list_epidemic_data[[i]]
    
    plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
         xlab = 'Time', ylab = ylabel,
         main = data_title,
         col = MODEL_COLOR,
         lwd = lwd, #3.5,
         cex.lab=cex+0.2, cex.axis=cex, cex.sub=cex-0.2,
         cex.main=cex+0.3, cex.main= main_font)   
  }
  
  #MAIN TITLE
  title(main_title, outer = TRUE, cex.main = main_font + 0.2)
  
  
  if(PDF){
    dev.off()
  }
}