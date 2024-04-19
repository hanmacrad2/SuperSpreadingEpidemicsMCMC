#********************************************
#1. PLOT_SIMULATION 
#********************************************
PLOT_SIMULATION_LIST <- function(list_epidemic_data, list_params,
                                 FLAGS_MODELS, FLAGS_PARAM, 
                                 MODEL_COLOR, RESULTS_FOLDER,
                                 #plot_width = 11.0, plot_height = 9.0,
                                 plot_width = 13.5, plot_height = 11.0,
                                    cex = 2.3, lwd = 2.45,
                                 main_font = 2.7, axis_font = 1.6,
                                    PDF = TRUE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAGS_PARAM)[which(unlist(FLAGS_PARAM))]
  ylabel = 'Infection count' #Daily
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_EPI_DATA_', param, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  #PLOT SETTINGS
  GET_MODEL_PLOT_SETTINGS(FLAGS_MODELS)
  
  #PLOT
  for(i in 1:length(list_epidemic_data)){
    
    #PARAMS
    #param = list_params[i]
    data_title = GET_MODEL_SIMULATION_TITLE(i, list_params, 
                                            FLAGS_MODELS, FLAGS_PARAM)
    #data_title = bquote(paste(.(model), " data. ", italic(R[0]),
    #                            " = ", .(param)))
    
    epidemic_data = list_epidemic_data[[i]]
    
    plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
         xlab = 'Time', ylab = ylabel,
         main = data_title,
         col = MODEL_COLOR,
         lwd = lwd, #3.5,
         cex.lab=cex+0.2, cex.axis=cex, cex.sub=cex-0.2,
         cex.main=cex+0.3, cex.main= main_font)   
  }
  
  
  if(PDF){
    dev.off()
  }
}

#GET_MODEL_SIMULATION_TITLE
GET_MODEL_SIMULATION_TITLE <- function(index, list_params, 
                                       FLAGS_MODELS, FLAGS_PARAM){
  
  #TITLE
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  if(FLAGS_MODELS$Baseline){
    param = list_params[i]
    data_title = bquote(paste(.(model), " data. ", italic(R[0]),
                                " = ", .(param)))
    #PLOT
    #par(mfrow = c(3,2))
    #par(mar = rep(4.5, 4)) #c(1.5, 5, 4, 1.5)) #bottom, left, top, right
    
  } else if (FLAGS_MODELS$SSE){
    
    #PARAMS
    #PLOT TITLE FIRST ROW
    if(index <= 3){
      r0_param = list_params$r0[index]
      k_param = list_params$k[index]
      
      if(FLAGS_PARAM$r0){
        data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                       " = ", bold(.(r0_param)),
                                       ", ", bold(k), " = ", bold(.(k_param))))) 
      } else if (FLAGS_PARAM$k){
        data_title = bquote(bold(paste(.(model), ": ", bold(k), " = ",
                                       bold(.(k_param)), ', ',
                                       bold(R[bold(0)]),
                                       " = ", bold(.(r0_param)))))
      }
      

    } else {
      data_title = ''
    }
    
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
  
  return(data_title)
}

#PLOT SETTINGS
GET_MODEL_PLOT_SETTINGS <- function(FLAGS_MODELS){
  
  if(FLAGS_MODELS$Baseline){
    
    #PLOT
    par(mfrow = c(3,2))
    par(mar = c(5, 5, 5, 5)) #rep(5, 4)) #c(1.5, 5, 4, 1.5)) #bottom, left, top, right
    
  } else if (FLAGS_MODELS$SSE){
    
    #PLOT
    par(mfrow = c(3,3))
    par(mar = c(4.9, 4.9, 5, 5)) 
   # par(mar = rep(4.8, 4))
    
  } else if (FLAGS_MODELS$SSI){
    
    #PLOT
    par(mfrow = c(3,3))
    
  } else if (FLAGS_MODELS$SSEB){
    
    #PLOT
    par(mfrow = c(3,3))
    
    
  } else if (FLAGS_MODELS$SSIB){
    
    #PLOT
    par(mfrow = c(3,3))
  }
}

#********************************************
#1. PLOT_SIMULATION 
#********************************************
PLOT_SIMULATION <- function(epidemic_data, FLAGS_MODELS, RESULTS_FOLDER, MODEL_COLOR,
                            plot_width = 10.0, plot_height = 9.0,
                            cex = 1.5, main_font = 2.5, axis_font = 1.6,
                            PDF = TRUE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  ylabel = 'Infection count' #Daily
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_EPI_DATA_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  #par(mar = c(5, 5, 4, 1.5)) #bottom, left, top, right
  
  if(FLAGS_MODELS$Baseline){
    par(mfrow = c(3,2))
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
  for(i in 1:6){
    plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
         xlab = 'Time', ylab = ylabel,
         main = data_title,
         col = MODEL_COLOR,
         lwd = 3.5, #3.5,
         cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)   
  }
  
  
  if(PDF){
    dev.off()
  }
}

#********************************************
#1. PLOT_SIMULATION 
#********************************************
PLOT_SIMULATION_LIST_V0 <- function(list_epidemic_data, 
                                 list_params, FLAGS_MODELS, RESULTS_FOLDER, MODEL_COLOR,
                            plot_width = 11.0, plot_height = 9.0,
                            cex = 1.65, main_font = 2.5, axis_font = 1.6,
                            PDF = TRUE){
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  ylabel = 'Infection count' #Daily
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_EPI_DATA_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  if(FLAGS_MODELS$Baseline){
    data_title = bquote(paste(.(model), " simulated epidemic data. ", italic(R[0]),
                              " = 2.0"))
    #PLOT
    par(mfrow = c(3,2))
    par(mar = rep(4.5, 4)) #c(1.5, 5, 4, 1.5)) #bottom, left, top, right
    
  } else if (FLAGS_MODELS$SSE){
    
    #PLOT
    par(mfrow = c(3,3))
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
  for(i in 1:length(list_epidemic_data)){
    
    #PARAMS
    param = list_params[i]
    data_title = 
    data_title = bquote(paste(.(model), " data. ", italic(R[0]),
                              " = ", .(param)))
    
    epidemic_data = list_epidemic_data[[i]]
    
    plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
         xlab = 'Time', ylab = ylabel,
         main = data_title,
         col = MODEL_COLOR,
         lwd = 3.5, #3.5,
         cex.lab=cex+0.2, cex.axis=cex, cex.main=cex+0.3, cex.sub=cex)   
  }
  
  
  if(PDF){
    dev.off()
  }
}

