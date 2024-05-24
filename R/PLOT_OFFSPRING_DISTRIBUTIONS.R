#********************************************
#1. PLOT_OFF SPRING DISTRIBUTIONS
#********************************************
PLOT_OFFSPRING_DISTRIBUTIONS <- function(list_params, FLAGS_MODELS, 
                                 MODEL_COLOR, RESULTS_FOLDER,
                                 num_offspring, 
                                 plot_width = 7.5, plot_height = 5.5,
                                 cex = 1.3, main_font = 1.85, #3.2, 
                                 axis_font = 1.6, PDF = TRUE){
  
  #PLOT SETTINGS
  #GET_OFFSPRING_TITLE(list_params, FLAGS_MODELS)
  
  #PARAMS
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  print(model)
  main_title =  bquote(paste(.(model), " model - Offspring Distribution"))
  param = list_params[1]
  
  if(PDF){
    #plot_folder = paste0(RESULTS_FOLDER, '/', model, '/')
    plot_folder = paste0(RESULTS_FOLDER, '/plot_num_offspring_', num_offspring, '/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_Offspring_dist_', param, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  #PLOT
  if(FLAGS_MODELS$SSEB){
    par(mar = c(5, 5, 5, 1)) #bottom left top right
  } else if(FLAGS_MODELS$SSIB){
    par(mar = c(4.5, 4.5, 3, 2)) #c(5, 4.7, 4, 4.7))
  } else {
    par(mar = rep(5,4))
  }
  
  #par(mar = c(5.5, 4.7, 4, 4.7))
  #par(oma = c(1, 1, 1, 1)) #bottom, left, top, right
  
  #DATA
  list_offspring = GET_OFFSPRING_DATA(list_params, FLAGS_MODELS, num_offspring)
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

#****************
#GET OFFSPRING DISTRIBUTIONS
GET_OFFSPRING_DATA <- function(list_params, FLAGS_MODELS, num_offspring){
  
  #DATA
  x <- 0:num_offspring  # Possible number of offspring
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[1]
    Z = dpois(x, lambda = r0_param)  # Compute PMF
    print(sum(x*Z))
     
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    
    r0_param = list_params[1]
    k_param = list_params[2]
    prob = k_param/(k_param + r0_param)
    Z = dnbinom(x, size = k_param, prob = prob)  
    print(sum(x*Z))
    
  } else if (FLAGS_MODELS$SSEB){
    Z = GET_OFFSPRING_SSEB(x, list_params, num_offspring)
    print(sum(x*Z))
    
  } else if (FLAGS_MODELS$SSIB){
    Z = GET_OFFSPRING_SSIB(x, list_params)
    print(sum(x*Z))
  }
  
  print(sum(x*Z))
  
  return(list(x = x, Z = Z))
}


#SSEB
GET_OFFSPRING_SSEB <-function(x, list_params, num_offspring, 
                              n_samps = 1000000){
  
  x <- 0:num_offspring
  # Parameters
  r0_param = list_params[1]
  alpha_param = list_params[2]
  beta_param = list_params[3]
  
  # Simulate X1 from the Poisson distribution with parameter alpha*R0
  pn = rpois(n_samps, lambda = alpha_param*r0_param)
  #pn <- dpois(x, lambda = alpha_param*r0_param) #sample
  
  # Simulate X2 from the compound Poisson distribution
  inner_lambda <- r0_param*(1 - alpha_param)/beta_param
  ps =  rpois(n_samps, lambda = beta_param*rpois(n_samps, lambda = inner_lambda))
  
  #ps <- beta_param*dpois(x, lambda = dpois(1, lambda = inner_lambda))
  p = pn + ps 
  
  #Density
  Z = rep(NA, length = num_offspring + 1)
  for (i in x){
    
    Z[i+1] = sum(p==i)
  }
  
  Z = Z/sum(Z) #probability vector
  
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


#SSEB
GET_OFFSPRING_SSEB_DATA <-function(list_params = c(2.5, 0.5, 10), n_samps = 100000,
                                   num_offspring = 50){
  
  x <- 0:num_offspring
  
  # Parameters
  r0_param = list_params[1]
  alpha_param = list_params[2]
  beta_param = list_params[3]
  
  # Simulate X1 from the Poisson distribution with parameter alpha*R0
  pn = rpois(n_samps, lambda = alpha_param*r0_param)
  #pn <- dpois(x, lambda = alpha_param*r0_param) #sample
  
  # Simulate X2 from the compound Poisson distribution
  inner_lambda <- r0_param*(1 - alpha_param)/beta_param
  ps =  rpois(n_samps, lambda = beta_param*rpois(n_samps, lambda = inner_lambda))
  
  #ps <- beta_param*dpois(x, lambda = dpois(1, lambda = inner_lambda))
  p = pn + ps 
  
  #Density
  Z = rep(NA, length = num_offspring + 1)
  for (i in x){
    
    Z[i+1] = sum(p==i)
  }
  
  Z = Z/sum(Z) #probability vector
  
  # Output the simulated total number of offspring Z
  return(Z)
  
}

#**************************
#* LEGEND 

#LEGEND 
GET_LEGEND_OFFSPRING <- function(list_params, MODEL_COLOR, FLAGS_MODELS,
                                 legend_location = 'topright', #alpha = 0.2,
                                 cex = 1.5, inset = -0.02){
  
  #PARAM
  legend_caption = GET_OFFSPRING_LEGEND_CAPTION(list_params, FLAGS_MODELS)
  
  #INSET
  if (FLAGS_MODELS$SSEB){
    inset = -0.05
  }
  
  #Legend
  legend(legend_location,
         x.intersp = 0.05, #-0.05,#CONTROLS RELATIVE TO THE PLOT
         legend_caption,
         cex = cex,
         inset = inset,
         #inset = c(-inset), #CONTROLS RELATIVE TO THE MARGINS NOT PLOT!
         col = c(MODEL_COLOR),
         lwd = 3, #rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = 1, # rep(1, num_conds), #c(1, 1),
         text.font = 2.3, #1.45
         bty = "n")
  
  #BIOMDAL MODELS EXTRA TEXT
  if (FLAGS_MODELS$SSEB){
    mtext(expression("|" ~ R[0] ~ "= 2.5 " ~ alpha ~ "= 0.5" ~ beta ~ "= 10"),
          side = 3, line = -6, cex = 1.35, adj = 0.95)

  } else if (FLAGS_MODELS$SSIB){
    mtext(expression("|" ~ R[0] ~ "=" ~ 2.5 ~ ", a = 0.5, b = 10"),
          side = 3, line = -6, cex = 1.35, adj = 0.95)
    #mtext("| R[0] = 2.5, a = 0.5, b = 10", side = 3, line = -6, cex = 1.25, adj = 0.95)
  }
}

#**************
# LEGEND CAPTION
GET_OFFSPRING_LEGEND_CAPTION <- function(list_params, FLAGS_MODELS){
  
  if(FLAGS_MODELS$Baseline){
    r0_param = list_params[1]
    #legend_offspring = bquote(paste('Poisson(', R[0], '= )', .(r0_param))) #DIDN'T WORK!!
    legend_offspring = expression(paste(' Z ~ Poisson(R'[0], '= 2.5)'))
    
  } else if (FLAGS_MODELS$SSE || FLAGS_MODELS$SSI){
    legend_offspring = expression(paste(' Z ~ NegBin(k, ', frac('k', paste('R'[0], ' + k')), ' ) | R'[0], '= 2.5, k = 0.4'))
    #legend_offspring = expression(paste('  NegBin(k, k/R'[0], '+ k )| R'[0], '= 2.5, k = 0.1'))
    
  } else if (FLAGS_MODELS$SSEB) {
    legend_offspring = expression(atop(
      paste('  Z ~ Poisson(', alpha * R[0], ')'),
      paste('+ Poisson(', beta, '* Poisson(R'[0] (1 - alpha)/beta ), '))'))
    
  }  else if (FLAGS_MODELS$SSIB){
    legend_offspring = expression(atop(
      paste('  Z ~ (1 - p) * Poisson(', a * R[0], ' + ', paste('(1 - a)R'[0]), '/', b, ')'),
      paste(' + ', p, '* Poisson(', a * b * R[0], ' + ', paste('(1 - a)R'[0]), ')')))
    
    legend_offspring2 = expression(
      atop(
        paste("Z ~ (1 - p) * Poisson(", a * R[0], " + ", (1 - a) * R[0] / b, ")"),
        paste("+ ", p, " * Poisson(", a * b * R[0], " + ", (1 - a) * R[0], ")"),
        atop(
          paste("| ", R[0], " = 2.5, a = 0.5, b = 10")
        )
      )
    )

  }
  
  return(legend_offspring)
}

