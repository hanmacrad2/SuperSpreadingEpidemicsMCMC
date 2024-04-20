#********************************************
#1. PLOT_SIMULATION 
#********************************************
PLOT_SIMULATION_LIST <- function(list_epidemic_data, list_params,
                                 FLAGS_MODELS, FLAGS_PARAM, 
                                 MODEL_COLOR, RESULTS_FOLDER,
                                 #plot_width = 11.0, plot_height = 9.0,
                                 plot_width = 13.5, plot_height = 11.0,
                                    cex = 2.3, lwd = 2.45,
                                 main_font = 3.2, axis_font = 1.6,
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

#GET_MODEL_SIMULATION_TITLE
GET_MODEL_SIMULATION_TITLE_V0 <- function(index, list_params, 
                                       FLAGS_MODELS, FLAGS_PARAM){
  
  #TITLE
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[index]
    data_title = bquote(bold(paste(bold(R[bold(0)]),
                                   " = ", bold(.(r0_param)))))
    main_title =  bquote(paste(.(model), " model - simulated epidemic data"))

    
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    
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
    
    if(index <= 7){ #7; beta
    
    r0_param = list_params$r0[index]
    alpha_param = list_params$alpha[index]
    beta_param = list_params$beta[index]
    
    if(FLAGS_PARAM$r0){
      data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                     " = ", bold(.(r0_param)),
                                     ", ", bold(alpha), " = ", bold(.(alpha_param)),
                                     ", ", bold(beta), " = ", bold(.(beta_param))))) 
    } else if (FLAGS_PARAM$alpha){
      data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                     " = ", bold(.(r0_param)),
                                     ", ", bold(alpha), " = ", bold(.(alpha_param)),
                                     ", ", bold(beta), " = ", bold(.(beta_param))))) 
    } else if (FLAGS_PARAM$beta){
      data_title = bquote(bold(paste(bold(beta), " = ",
                                     bold(.(beta_param))))) #, ', ',
                                     # bold(R[bold(0)]),
                                     # " = ", bold(.(r0_param)),
                                     # ', ',
                                     # bold(alpha),
                                     # " = ", bold(.(alpha_param)))))
    }
    
    } else {
      data_title = ''
    }
    
  } else if (FLAGS_MODELS$SSIB){
    
    r0_param = list_params$r0[index]
    a_param = list_params$a[index]
    b_param = list_params$b[index]
    
    data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                   " = ", bold(.(r0_param)),
                                   ", ", bold(a), " = ", bold(.(a_param)),
                                   ", ", bold(b), " = ", bold(.(b_param)))))
  }
  
  list_titles = list(data_title= data_title, main_title = main_title)
  
  return(list_titles)
}

#PLOT SETTINGS
GET_MODEL_PLOT_SETTINGS <- function(FLAGS_MODELS){
  
  if(FLAGS_MODELS$Baseline){
    
    #PLOT
    par(mfrow = c(3,2))
    par(mar = c(5.5, 4.7, 4, 4.7)) #rep(4.5, 4)) #c(1.5, 5, 4, 1.5)) #bottom, left, top, right
    par(oma = c(1, 1, 6, 1))
    
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    
    #PLOT
    par(mfrow = c(3,3))
    par(mar = c(5.5, 5, 4, 4.7)) #rep(4.5, 4)) #c(1.5, 5, 4, 1.5)) #bottom, left, top, right
    par(oma = c(1, 1, 6, 1))
    #par(mar = c(4.9, 4.9, 5, 5)) 
   # par(mar = rep(4.8, 4))
    
  } else if (FLAGS_MODELS$SSI){
    
    #PLOT
    #par(mfrow = c(3,3))
    
  } else if (FLAGS_MODELS$SSEB){
    
    #PLOT
    par(mfrow = c(3,2))
    par(mar = c(5, 5, 5, 5))
    
  } else if (FLAGS_MODELS$SSIB){
    
    #PLOT
    par(mfrow = c(3,2))
    par(mar = c(5, 5, 5, 5))
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
#* PLOT SSIB SIMULATION LIST
PLOT_SSIB_SIMULATION_SS_LIST <- function(list_ssib_data_total, RESULTS_FOLDER,
                                    list_ssib_params = list(r0 = 2, a = 0.5, b = 10),
                                    model = 'SSIB', 
                                    col_ss = 'aquamarine', col_ns = 'orange',
                                    plot_width = 14.2, plot_height = 11.0,
                                    cex = 2.3, lwd_data = 2.4,
                                    main_font = 2.5, axis_font = 1.6,
                                    PDF = TRUE){
  
  #PLOT
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_EPI_DATA_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  par(mfrow = c(2,2))
  par(mar = c(4.25, 4.7, 4, 4.2))
  #par(mar = c(5, 5, 4, 1.5)) #bottom, left, top, right
  
  #PLOT
  for(i in 1:length(list_ssib_data_total)){
    
    list_ssib_data = list_ssib_data_total[[i]]
    
    epidemic_data = list_ssib_data$total_infections
    non_ss = list_ssib_data$non_ss
    ss = list_ssib_data$ss
    ylabel = 'Infection count' 
    
    #PLOT
    r0_param = list_ssib_params$r0[i]
    a_param = list_ssib_params$a[i]
    b_param = list_ssib_params$b[i]
    
    data_title = bquote(bold(paste(.(model), ": ", bold(R[0]),
                                   " = ", .(r0_param), ", ", bold(a), " = ", .(a_param),
                                   ", ", bold(b), ' = ', .(b_param))))
    
    #PLOT
    plot(seq_along(epidemic_data), epidemic_data,  type = 'l',
         xlab = 'Time', ylab = ylabel,
         main = data_title,
         lwd = lwd_data, 
         cex.lab=cex, cex.axis=cex, 
        cex.sub=cex, cex.main= main_font)  #cex.main=cex+0.2, 
    
    #NS 
    lines(seq_along(non_ss), non_ss, col = col_ns, lwd = lwd_data)
    
    #SS
    lines(seq_along(ss), ss, col = col_ss, lwd = lwd_data)
    
    legend('topleft', c('Epidemic data in total', 'Non super-spreaders', 'Super-spreaders'),
           col = c('black', col_ns, col_ss), 
           lwd = c(2, 2, 2),
           cex = 1.6,
           lty = c(1,1,1))
    
  }
  
  
  if(PDF){
    dev.off()
  }
}

#********************************************
#GET_MODEL_SIMULATION_TITLE
GET_MODEL_SIMULATION_TITLE <- function(index, list_params, 
                                       FLAGS_MODELS, FLAGS_PARAM){
  
  #TITLE
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[index]
    data_title = bquote(bold(paste(bold(R[bold(0)]),
                                   " = ", bold(.(r0_param)))))
    main_title =  bquote(paste(.(model), " model - simulated epidemic data"))
    
    
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    
    #TITLE
    r0_param = list_params$r0[index]
    k_param = list_params$k[1]
    main_title = bquote(paste(.(model), " model epidemic data - varying ",
                              R[0], ', constant k = ', .(k_param)))
    if(index <= 3){
      r0_param = list_params$r0[index]
      k_param = list_params$k[index]
      
      if(FLAGS_PARAM$r0){
        
        data_title = bquote(bold(paste(bold(R[bold(0)]),
                                       " = ", bold(.(r0_param)))))
                                      # ", ", bold(k), " = ", bold(.(k_param)))))
        #data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
        #                               " = ", bold(.(r0_param)),
        #                               ", ", bold(k), " = ", bold(.(k_param)))))
        
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
    
    if(index <= 7){ #7; beta
      
      r0_param = list_params$r0[index]
      alpha_param = list_params$alpha[index]
      beta_param = list_params$beta[index]
      
      if(FLAGS_PARAM$r0){
        data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                       " = ", bold(.(r0_param)),
                                       ", ", bold(alpha), " = ", bold(.(alpha_param)),
                                       ", ", bold(beta), " = ", bold(.(beta_param))))) 
      } else if (FLAGS_PARAM$alpha){
        data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                       " = ", bold(.(r0_param)),
                                       ", ", bold(alpha), " = ", bold(.(alpha_param)),
                                       ", ", bold(beta), " = ", bold(.(beta_param))))) 
      } else if (FLAGS_PARAM$beta){
        data_title = bquote(bold(paste(bold(beta), " = ",
                                       bold(.(beta_param))))) #, ', ',
      }
      
    } else {
      data_title = ''
    }
    
  } else if (FLAGS_MODELS$SSIB){
    
    r0_param = list_params$r0[index]
    a_param = list_params$a[index]
    b_param = list_params$b[index]
    
    data_title = bquote(bold(paste(.(model), ": ", bold(R[bold(0)]),
                                   " = ", bold(.(r0_param)),
                                   ", ", bold(a), " = ", bold(.(a_param)),
                                   ", ", bold(b), " = ", bold(.(b_param)))))
  }
  
  list_titles = list(data_title= data_title, main_title = main_title)
  
  return(list_titles)
}