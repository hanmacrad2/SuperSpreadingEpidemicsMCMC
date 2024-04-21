#********************************************
#1. PLOT_OFF SPRING DISTRIBUTIONS
#********************************************
PLOT_OFFSPRING_DISTRIBUTIONS <- function(list_params, FLAGS_MODELS, 
                                 MODEL_COLOR, RESULTS_FOLDER,
                                 num_offspring = 25, r0 = 2, 
                                 plot_width = 7, plot_height = 5.0,
                                 cex = 1.3, main_font = 2.0, #3.2, 
                                 axis_font = 1.6, PDF = TRUE){
  
  #PLOT SETTINGS
  main_title = GET_OFFSPRING_TITLE(list_params, FLAGS_MODELS)
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = list_params[1]
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_Offspring_dist_', param, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  #DATA
  list_offspring = GET_OFFSPRING_DATA(list_params, FLAGS_MODELS)
  x = list_offspring$x
  pmf_offspring = list_offspring$pmf
  
  barplot(pmf_offspring, names.arg = x, 
          col = MODEL_COLOR, #"skyblue",
          main = main_title,
          xlab = "Number of Offspring",
          ylab = "Probability",
          #lwd = lwd, #3.5,
          cex.lab=cex+0.2, cex.axis=cex, cex.sub=cex-0.2,
          cex.main=cex+0.3, cex.main= main_font) 
  
  #GET_MODEL_PLOT_SETTINGS(FLAGS_MODELS, FLAGS_PARAM)

  #MAIN TITLE
  #title(main_title, outer = TRUE, cex.main = main_font + 0.2)
  
  if(PDF){
    dev.off()
  }
}

GET_OFFSPRING_TITLE <- function(list_params, FLAGS_MODELS){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[1]
    main_title =  bquote(paste(.(model), " model - Offspring Distriubtion. Poisson(", R[0], ' = ', .(r0_param), ')'))
  
    } else if (FLAGS_MODELS$SSE){
    
    }
  
  return(main_title)
}


GET_OFFSPRING_DATA <- function(list_params, FLAGS_MODELS, num_offspring = 25){
  
  #DATA
  x <- 0:num_offspring  # Possible number of offspring
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[1]
    pmf = dpois(x, lambda = r0_param)  # Compute PMF
     
  } else if (FLAGS_MODELS$SSE){
    
  }
  
  return(list(x = x, pmf = pmf))
}







# #PLOT
# plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
#      xlab = 'Time', ylab = ylabel,
#      main = data_title,
#      col = MODEL_COLOR,
#      lwd = lwd, #3.5,
#      cex.lab=cex+0.2, cex.axis=cex, cex.sub=cex-0.2,
#      cex.main=cex+0.3, cex.main= main_font)  
# 
# #PLOTS
# R0 <- 2  
# cex = 1.7
# # Generate data
# x <- 0:25  # Possible number of offspring
# pmf <- dpois(x, lambda = R0)  # Compute PMF
# 
# # Plot
# model = 'Baseline'
# main_title =  bquote(paste(.(model), " model - Offspring Distriubtion. Poisson(", R[0], ' = 2)'))
# 
# # Generate data
# x <- 0:25  # Possible number of offspring
# pmf <- dpois(x, lambda = R0)  # Compute PMF
# 
# barplot(pmf, names.arg = x, 
#         col = MODEL_COLORS[1], #"skyblue",
#         main = main_title,
#         xlab = "Number of Offspring",
#         ylab = "Probability",
#         lwd = lwd, #3.5,
#         cex.lab=cex+0.2, cex.axis=cex, cex.sub=cex-0.2,
#         cex.main=cex+0.3, cex.main= main_font) 
