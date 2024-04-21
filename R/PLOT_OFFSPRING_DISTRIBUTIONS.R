#********************************************
#1. PLOT_OFF SPRING DISTRIBUTIONS
#********************************************
PLOT_OFFSPRING_DISTRIBUTIONS <- function(list_params, FLAGS_MODELS, 
                                 MODEL_COLOR, RESULTS_FOLDER,
                                 num_offspring = 25, r0 = 2, 
                                 plot_width = 7.5, plot_height = 5.5,
                                 cex = 1.3, main_font = 1.85, #3.2, 
                                 axis_font = 1.6, PDF = TRUE){
  
  #PLOT SETTINGS
  #GET_OFFSPRING_TITLE(list_params, FLAGS_MODELS)
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  main_title =  bquote(paste(.(model), " model - Offspring Distriubtion"))
  param = list_params[1]
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_Offspring_dist_', param, '_', time_stamp, '.pdf')  
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
  
  #DATA
  list_offspring = GET_OFFSPRING_DATA(list_params, FLAGS_MODELS)
  x = list_offspring$x
  Z = list_offspring$Z
  
  barplot(Z, names.arg = x, 
          col = MODEL_COLOR, #"skyblue",
          main = main_title,
          xlab = "Number of Offspring (Z)",
          ylab = "Probability",
          #lwd = lwd, #3.5,
          cex.lab=cex + 0.25, cex.axis=cex + 0.1, 
          cex.main=cex+0.3, cex.main= main_font) 

  GET_LEGEND_OFFSPRING(list_params, MODEL_COLOR, FLAGS_MODELS)

  #GET_MODEL_PLOT_SETTINGS(FLAGS_MODELS, FLAGS_PARAM)

  #MAIN TITLE
  #title(main_title, outer = TRUE, cex.main = main_font + 0.2)
  
  if(PDF){
    dev.off()
  }
}

#DATA
GET_OFFSPRING_DATA <- function(list_params, FLAGS_MODELS, num_offspring = 25){
  
  #DATA
  x <- 0:num_offspring  # Possible number of offspring
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[1]
    Z = dpois(x, lambda = r0_param)  # Compute PMF
     
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    
    r0_param = list_params[1]
    k_param = list_params[2]
    prob = k_param/(k_param + r0_param)
    Z = dnbinom(x, size = k_param, prob = prob)  
    
  } else if (FLAGS_MODELS$SSEB){
    Z = GET_OFFSPRING_SSEB(x, list_params)
    
  } else if (FLAGS_MODELS$SSIB){
    Z = GET_OFFSPRING_SSIB(x, list_params)
  }
  
  return(list(x = x, Z = Z))
}


#SSEB
GET_OFFSPRING_SSEB <-function(x, list_params){
  
  # Parameters
  r0_param = list_params[1]
  alpha_param = list_params[2]
  beta_param = list_params[3]
  
  # Simulate X1 from the Poisson distribution with parameter alpha*R0
  pois1_nsei <- dpois(x, lambda = alpha_param*r0_param)
  
  # Simulate X2 from the compound Poisson distribution
  inner_lambda <- r0_param*(1 - alpha_param)/beta_param
  pois2_ssei <- beta_param*dpois(x, lambda = dpois(1, lambda = inner_lambda))
  
  # Calculate the total sum Z
  Z <- pois1_nsei + pois2_ssei
  
  # Output the simulated total number of offspring Z
  return(Z)
  
}

GET_OFFSPRING_SSIB <- function(x, list_params){
  
  # Parameters
  # Parameters
  r0 = list_params[1]
  a = list_params[2]
  b = list_params[3]
  p = (1 - a)/(1 - a + a*b)
  
  # Calculate the mean of the first Poisson distribution
  lambda1 <- a*r0 + ((1 - a)*r0)/b
  # Calculate the mean of the second Poisson distribution
  lambda2 <- a*b*r0 + (1 - a)*r0
  
  # Simulate from the first Poisson distribution
  Z1 <- dpois(x, lambda = lambda1)
  
  # Simulate from the second Poisson distribution
  Z2 <- dpois(x, lambda = lambda2)
  
  # Weighted sum of the two distributions
  Z <- (1-p)*Z1 + p*Z2
  
  # Output the simulated total number of offspring Z
  return(Z)
  
}

#CAPTION
GET_OFFSPRING_LEGEND_CAPTION <- function(list_params, FLAGS_MODELS){
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[1]
    #legend_offspring = bquote(paste('Poisson(', R[0], '= )', .(r0_param))) #DIDN'T WORK!!
    legend_offspring = expression(paste(' Z ~ Poisson(R'[0], '= 2.5)'))
    
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    legend_offspring = expression(paste(' Z ~ NegBin(k, ', frac('k', paste('R'[0], ' + k')), ' ) | R'[0], '= 2.5, k = 0.1'))
    #legend_offspring = expression(paste('  NegBin(k, k/R'[0], '+ k )| R'[0], '= 2.5, k = 0.1'))
  
  } else if (FLAGS_MODELS$SSEB) {
    legend_offspring = expression(atop(
      paste('  Z ~ Poisson(', alpha * R[0]),
      paste(' + ', beta, '* Poisson()', paste('(1 - a)R'[0]), ')')))
    
  }  else if (FLAGS_MODELS$SSIB){
    #legend_offspring = expression(paste('  Z ~ 1-p*Poisson(aR'[0], ' +  ', frac(paste('(1 - a)R'[0]), b), ' )'))
                                       # , ', frac('k', paste('R'[0], ' + k')), ' ) | R'[0], '= 2.5, k = 0.1')) 
    
      legend_offspring = expression(atop(
      paste('  Z ~ (1 - p) * Poisson(', a * R[0], ' + ', paste('(1 - a)R'[0]), '/', b, ')'),
      paste(' + ', p, '* Poisson(', a * b * R[0], ' + ', paste('(1 - a)R'[0]), ')')))
      #paste('| ', R[0], ' = 2.5, a = 0.5, b = 10')))
    
    legend_offspring_wrong = expression(paste(
      ' Z ~ (1 - p) * Poisson(', a * R[0], ' + ', paste('(1 - a)R'[0]), '/', b, ') \n',
      p, '* Poisson(', a * b * R[0], ' + ', paste('(1 - a)R'[0]), ') \n',
      '| ', R[0], ' = 2.5, a = 0.5, b = 10'
    ))
    
    
    
    }
  
  return(legend_offspring)
}

#LEGEND 
GET_LEGEND_OFFSPRING <- function(list_params, MODEL_COLOR, FLAGS_MODELS,
                            legend_location = 'topright', #alpha = 0.2,
                            cex = 1.5, inset = 0.1){
  
  #PARAM
  legend_caption = GET_OFFSPRING_LEGEND_CAPTION(list_params, FLAGS_MODELS)
  
  #Legend
  legend(legend_location,
         x.intersp = 0.05, #-0.05,#CONTROLS RELATIVE TO THE PLOT
         legend_caption,
         cex = cex,
         inset = -0.02,
         #inset = c(-inset), #CONTROLS RELATIVE TO THE MARGINS NOT PLOT!
         col = c(MODEL_COLOR),
         lwd = 3, #rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = 1, # rep(1, num_conds), #c(1, 1),
         text.font = 2.3, #1.45
         bty = "n")
}






# GET_OFFSPRING_TITLE <- function(list_params, FLAGS_MODELS){
#   
#   #PARAMS
#   model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
#   
#   if(FLAGS_MODELS$Baseline){
#     r0_param = list_params[1]
#     main_title =  bquote(paste(.(model), " model - Offspring Distriubtion")) #. Poisson(", R[0], ' = ', .(r0_param), ')'))
#     #main_title =  bquote(paste(.(model), " model - Offspring Distriubtion. Poisson(", R[0], ' = ', .(r0_param), ')'))
#     
#   } else if (FLAGS_MODELS$SSE){
#     
#   }
#   
#   return(main_title)
# }

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
