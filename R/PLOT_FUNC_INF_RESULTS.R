#PLOT INFERENCE RESULTS
PLOT_CI_PARAMS <- function(df_results, COMP_FOLDER, fig_num = '1',
                           num_days = 50, cex = 1.0, num_conds = 5,
                             PDF = FALSE, FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                               alpha = FALSE, gamma = FALSE),
                             FLAG_FILTER = list(end_day = FALSE,
                                                tot_infs = FALSE), 
                           FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                             SSEB = FALSE, SSIB = FALSE)){

  #Params
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  true_total = unlist(df_results[paste0('true_', true_param)])
  filter_param = names(FLAG_FILTER)[which(unlist(FLAG_FILTER))]
  model =  names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  #paste0(model, ' Model ', true_param, ' Inference')
  
  #Plot
  if(PDF){
    pdf_file = paste0('Fig_', model, '_param_', true_param, '_', fig_num, '.pdf') 
    pdf(paste0(COMP_FOLDER, '/plots/', pdf_file), width = 13.0, height = 8.0)
  }
  par(mar=c(4.9, 4.1, 3.0, 15.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  #Colours
  viridis_palette <- viridis(10)
  selected_colors <- viridis_palette[c(9, 10, 6, 8, 5)]
  
  # Define different filters and their corresponding labels
  filter_list <- list(
    df2 = df_results[df_results[[filter_param]] == 2, ],
    df3 = df_results[df_results[[filter_param]] == 3, ],
    df4 = df_results[(df_results[[filter_param]] == 4) | (df_results[[filter_param]] == 5), ],
    df5 = df_results[(df_results[[filter_param]] > 5) & (df_results[[filter_param]] <= 10), ],
    df10 = df_results[df_results[[filter_param]] > 10, ]
  )
  
  # Create a list of legends for each subset
  if(FLAG_FILTER$tot_infs){
    legend_list <- c(paste0(true_param, " True "),
                     paste0(true_param, " posterior mean, infs: [2]"),
                     paste0(true_param, " post. mean, infs: [3]"),
                     paste0(true_param, " post. mean, infs: [4, 5]"),
                     paste0(true_param, " post. mean, infs: [6, 10]"),
                     paste0(true_param, " post. mean, infs > 10"))
    }
  y_lim = c(0, max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
  x_lim = c(min(df_results[paste0('true_', true_param)]), max(df_results[paste0('true_', true_param)]))
  
  # LABELS
  if (FLAG_PARAM$r0) {
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
    ylab <- expression(paste('Posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX = bquote(bold(paste(.(model) ~ "Model - ", italic(R[0]), " Inference. Epidemic length: " ~ .(num_days) ~ ' days')))
      
  } else {
    xlab <- paste0('True ', true_param)
    ylab <- paste0('Posterior mean of ', true_param)
    titleX = bquote(bold(paste(.(model) ~ "Model - " ~ .(true_param) ~
                                 "Inference. Epidemic length: " ~ .(num_days) ~ ' days')))
  }
  
  # Create the plot for each subset
  for (i in 1:length(filter_list)) {
    
    df_subset <- filter_list[[i]]
    print(paste0('n row:', nrow(df_subset)))
    true_subset <- unlist(df_subset[paste0('true_', true_param)])
    mean_subset <- unlist(df_subset[paste0('mean_', true_param)])
    upper_ci_subset <- unlist(df_subset[paste0('upper_ci_', true_param)])
    lower_ci_subset <- unlist(df_subset[paste0('lower_ci_', true_param)])
    colour <- selected_colors[i]
      
    if (i == 1) {
      
      plot(true_subset, mean_subset, type = "p",
           main = "",
           xlab = xlab, ylab = ylab,
           xlim = x_lim,
           ylim = y_lim,
           col = colour, pch = 16, 
           #text.font = 4.0,
           cex = cex)
      title(main = list(titleX, cex = 1.2, font = 2.0))
      
    } else {
      
      points(true_subset, mean_subset, type = "p",
           xlab = xlab, ylab = ylab,
           ylim = c(0, max(true_subset, upper_ci_subset)),
           col = colour, pch = 16, cex = cex)
    }
    
    # Add error bars using segments()
    segments(true_subset, upper_ci_subset, true_subset,
             lower_ci_subset, lwd = 0.5, col = colour)
    
    # Points
    points(true_subset, mean_subset, type = "p",
           col = colour, pch = 16)
    
  }
  
  # TRUE line
  lines(true_total, true_total, col = 'black', lwd = 3)
  
  # Legend
  legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         inset=c(-0.32,0),
         col = c('black', selected_colors),
         lwd = rep(1.8, num_conds),
         lty = rep(1, num_conds), #c(1, 1),
         pch = c(NA, 19, 19, 19, 19, 19),
         text.font = 1.0,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}

#*********************
#* PLOT 2
PLOT_CI_PARAMS_II <- function(df_results, COMP_FOLDER, fig_num = '1',
                           cex = 1.0, num_conds = 5,
                           PDF = FALSE, 
                           FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                             alpha = FALSE, gamma = FALSE),
                           FLAG_FILTER = list(end_day = FALSE,
                                              tot_infs = FALSE), 
                           FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                             SSEB = FALSE, SSIB = FALSE)){
  
  #Params
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  true_total = unlist(df_results[paste0('true_', true_param)])
  filter_param = names(FLAG_FILTER)[which(unlist(FLAG_FILTER))]
  model =  names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  #paste0(model, ' Model ', true_param, ' Inference')
  
  #Plot
  if(PDF){
    pdf_file = paste0('Fig_', model, '_param_', true_param, '_', fig_num, '.pdf') 
    pdf(paste0(COMP_FOLDER, '/plots/', pdf_file), width = 13.0, height = 8.0)
  }
  par(mar=c(4.9, 4.1, 3.0, 15.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  #Colours
  viridis_palette <- viridis(10)
  #selected_colors <- viridis_palette[c(9, 10, 6, 8, 5)]
  selected_colors <- viridis_palette[c(9, 10, 8, 5)]
  
  # Define different filters and their corresponding labels
  #browser()
  filter_list <- list(
    df5 = df_results[(df_results[[filter_param]] > 4) & (df_results[[filter_param]] <= 10), ],
    df10 = df_results[(df_results[[filter_param]] > 10) & (df_results[[filter_param]] <= 100), ],
    df100 = df_results[(df_results[[filter_param]] > 100) & (df_results[[filter_param]] <= 500), ],
    df500 = df_results[df_results[[filter_param]] > 500, ]
  )
  
  # Create a list of legends for each subset
  if(FLAG_FILTER$tot_infs){
    legend_list <- c(paste0(true_param, " True "),
                     paste0(true_param, " post. mean, infs: [5, 10]"),
                     paste0(true_param, " post. mean, infs: [10, 100]"),
                     paste0(true_param, " post. mean, infs: [100, 500]"),
                     paste0(true_param, " post. mean, infs > 500"))
  }
  y_lim = c(0, max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
  x_lim = c(min(df_results[paste0('true_', true_param)]), max(df_results[paste0('true_', true_param)]))
  
  # LABELS
  if (FLAG_PARAM$r0) {
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
    ylab <- expression(paste('Posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX = bquote(bold(paste(.(model) ~ "Model - ", italic(R[0]), " Inference")))
    
  } else {
    xlab <- paste0('True ', true_param)
    ylab <- paste0('Posterior mean of ', true_param)
    titleX = bquote(bold(paste(.(model) ~ "Model - " ~ .(true_param) ~ "Inference")))
  }
  
  # Create the plot for each subset
  for (i in 1:length(filter_list)) {
    
    df_subset <- filter_list[[i]]
    print(paste0('n row:', nrow(df_subset)))
    true_subset <- unlist(df_subset[paste0('true_', true_param)])
    mean_subset <- unlist(df_subset[paste0('mean_', true_param)])
    upper_ci_subset <- unlist(df_subset[paste0('upper_ci_', true_param)])
    lower_ci_subset <- unlist(df_subset[paste0('lower_ci_', true_param)])
    colour <- selected_colors[i]
    
    if (i == 1) {
      
      plot(true_subset, mean_subset, type = "p",
           main = "", xlab = xlab, ylab = ylab,
           xlim = x_lim, ylim = y_lim,
           col = colour, pch = 16, cex = cex)#text.font = 4.0,
      title(main = list(titleX, cex = 1.2, font = 2.0))
      
    } else {
      
      points(true_subset, mean_subset, type = "p",
             xlab = xlab, ylab = ylab,
             ylim = c(0, max(true_subset, upper_ci_subset)),
             col = colour, pch = 16, cex = cex)
    }
    
    # Add error bars using segments()
    segments(true_subset, upper_ci_subset, true_subset,
             lower_ci_subset, lwd = 0.5, col = colour)
    
    # Points
    points(true_subset, mean_subset, type = "p",
           col = colour, pch = 16)
    
  }
  
  # TRUE line
  lines(true_total, true_total, col = 'black', lwd = 3)
  
  # Legend
  legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         inset=c(-0.32,0),
         col = c('black', selected_colors),
         lwd = rep(1.8, num_conds),
         lty = rep(1, num_conds), #c(1, 1),
         pch = c(NA, 19, 19, 19, 19, 19),
         text.font = 1.0,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}

#PLOT INFERENCE RESULTS
PLOT_CI_PARAMS_FIXED <- function(df_results, COMP_FOLDER, fig_num = '1',
                           num_days = 50, cex = 1.0, num_conds = 5,
                           PDF = FALSE, 
                           FIXED_PARAM = list(r0 = FALSE, k = FALSE,
                                             alpha = FALSE, gamma = FALSE),
                           VAR_PARAM = list(r0 = FALSE, k = FALSE,
                                              alpha = FALSE, gamma = FALSE),
                           FLAG_FILTER = list(end_day = FALSE,
                                              tot_infs = TRUE), 
                           FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                             SSEB = FALSE, SSIB = FALSE)){

  #Params
  var_param =  names(VAR_PARAM)[which(unlist(VAR_PARAM))]
  fixed_param = names(FIXED_PARAM)[which(unlist(FIXED_PARAM))]
  fixed_total = unlist(df_results[paste0('true_', fixed_param)])
  var_total = unlist(df_results[paste0('true_', var_param)])
  filter_param = names(FLAG_FILTER)[which(unlist(FLAG_FILTER))]
  model =  names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  #paste0(model, ' Model ', true_param, ' Inference')
  
  #Plot
  if(PDF){
    pdf_file = paste0('Fig_', model, '_fixed_', fixed_param, '_', fig_num, '.pdf')  
    pdf(paste0(COMP_FOLDER, '/plots/', pdf_file), width = 13.0, height = 8.0)
  }
  
  par(mar=c(4.9, 4.1, 3.0, 15.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  #Colours
  viridis_palette <- viridis(10)
  selected_colors <- viridis_palette[c(9, 10, 6, 8, 5)]
  
  # Define different filters and their corresponding labels
  filter_list <- list(
    df2 = df_results[df_results[[filter_param]] == 2, ],
    df3 = df_results[df_results[[filter_param]] == 3, ],
    df4 = df_results[(df_results[[filter_param]] == 4) | (df_results[[filter_param]] == 5), ],
    df5 = df_results[(df_results[[filter_param]] > 5) & (df_results[[filter_param]] <= 10), ],
    df10 = df_results[df_results[[filter_param]] > 10, ]
  )
  
  # Create a list of legends for each subset
  if(FLAG_FILTER$tot_infs){
    print('LEGEND')
    legend_list <- c(paste0(fixed_param, " True "),
                     paste0(fixed_param, " posterior mean, infs: [2]"),
                     paste0(fixed_param, " post. mean, infs: [3]"),
                     paste0(fixed_param, " post. mean, infs: [4, 5]"),
                     paste0(fixed_param, " post. mean, infs: [6, 10]"),
                     paste0(fixed_param, " post. mean, infs > 10"))
    print(legend_list)
  }
  y_lim = c(0, max(fixed_total, max(df_results[paste0('upper_ci_', fixed_param)]))) 
  x_lim = c(min(df_results[paste0('true_', var_param)]), max(df_results[paste0('true_', var_param)]))
  
  # LABELS
  if (FIXED_PARAM$r0){
    ylab <- expression(paste('Posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX = bquote(bold(paste(.(model) ~ "Model - ", italic(R[0]), " Inference. Epidemic length: " ~ .(num_days) ~ ' days')))
  } else {
    ylab <- paste0('Posterior mean of ', fixed_param)
    titleX = bquote(bold(paste(.(model) ~ "Model - " ~ .(fixed_param) ~  " Inference. Epidemic length: " ~ .(num_days) ~ ' days')))
  }
  
  #x label
  if (VAR_PARAM$r0) { #Var param should be x axis 
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
  } else {
    xlab <- paste0('True ', var_param)
  }
  
  # Create the plot for each subset
  for (i in 1:length(filter_list)) {
    
    df_subset <- filter_list[[i]]
    print(paste0('n row:', nrow(df_subset)))
    var_subset <- unlist(df_subset[paste0('true_', var_param)])
    mean_subset <- unlist(df_subset[paste0('mean_', fixed_param)])
    upper_ci_subset <- unlist(df_subset[paste0('upper_ci_', fixed_param)])
    lower_ci_subset <- unlist(df_subset[paste0('lower_ci_', fixed_param)])
    colour <- selected_colors[i]
    
    if (i == 1) {
      
      plot(var_subset, mean_subset, type = "p",
           main = "",
           xlab = xlab, ylab = ylab,
           xlim = x_lim,
           ylim = y_lim,
           col = colour, pch = 16, 
           #text.font = 4.0,
           cex = cex)
      title(main = list(titleX, cex = 1.2, font = 2.0))
      
    } else {
      
      points(var_subset, mean_subset, type = "p",
             xlab = xlab, ylab = ylab,
             ylim = c(0, max(var_subset, upper_ci_subset)),
             col = colour, pch = 16, cex = cex)
    }
    
    #Error bars 
    segments(var_subset, upper_ci_subset, var_subset,
             lower_ci_subset, lwd = 0.5, col = colour)
    # Points
    points(var_subset, mean_subset, type = "p",
           col = colour, pch = 16)
  }
  # TRUE line
  lines(var_total, fixed_total, col = 'black', lwd = 3)
  
  # Legend
  legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         inset=c(-0.32,0),
         col = c('black', selected_colors),
         lwd = rep(1.8, num_conds),
         lty = rep(1, num_conds), #c(1, 1),
         pch = c(NA, 19, 19, 19, 19, 19),
         text.font = 1.0,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}

PLOT_CI_PARAMS_FIXED_II <- function(df_results, COMP_FOLDER, fig_num = '1',
                                 cex = 1.0, num_conds = 5,
                                 PDF = FALSE, 
                                 FIXED_PARAM = list(r0 = FALSE, k = FALSE,
                                                    alpha = FALSE, gamma = FALSE, c = FALSE),
                                 VAR_PARAM = list(r0 = FALSE, k = FALSE,
                                                  alpha = FALSE, gamma = FALSE, c = FALSE),
                                 FLAG_FILTER = list(end_day = FALSE,
                                                    tot_infs = TRUE), 
                                 FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                                   SSEB = FALSE, SSIB = FALSE)){
  
  #Params
  var_param =  names(VAR_PARAM)[which(unlist(VAR_PARAM))]
  fixed_param = names(FIXED_PARAM)[which(unlist(FIXED_PARAM))]
  fixed_total = unlist(df_results[paste0('true_', fixed_param)])
  var_total = unlist(df_results[paste0('true_', var_param)])
  filter_param = names(FLAG_FILTER)[which(unlist(FLAG_FILTER))]
  model =  names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  #paste0(model, ' Model ', true_param, ' Inference')
  
  #Plot
  if(PDF){
    pdf_file = paste0('Fig_', model, '_fixed_', fixed_param, '_', fig_num, '.pdf')  
    pdf(paste0(COMP_FOLDER, '/plots/', pdf_file), width = 13.5, height = 8.0)
  }
  
  par(mar=c(4.9, 4.1, 3.0, 15.5), xpd=TRUE) #Margins; bottom, left, top, right
  
  #Colours
  viridis_palette <- viridis(10)
  selected_colors <- viridis_palette[c(9, 10, 8, 5)]
  #selected_colors <- viridis_palette[c(9, 10, 6, 8, 5)]
  
  # Define different filters and their corresponding labels
  filter_list <- list(
    df5 = df_results[(df_results[[filter_param]] > 4) & (df_results[[filter_param]] <= 10), ],
    df10 = df_results[(df_results[[filter_param]] > 10) & (df_results[[filter_param]] <= 100), ],
    df100 = df_results[(df_results[[filter_param]] > 100) & (df_results[[filter_param]] <= 500), ],
    df500 = df_results[df_results[[filter_param]] > 500, ]
  )
  
  # Create a list of legends for each subset
  if(FLAG_FILTER$tot_infs){
    legend_list <- c(paste0(fixed_param, " True "),
                     paste0(fixed_param, " post. mean, infs: [5, 10]"),
                     paste0(fixed_param, " post. mean, infs: [10, 100]"),
                     paste0(fixed_param, " post. mean, infs: [100, 500]"),
                     paste0(fixed_param, " post. mean, infs > 500"))
  }

  y_lim = c(0, max(fixed_total, max(df_results[paste0('upper_ci_', fixed_param)]))) 
  x_lim = c(min(df_results[paste0('true_', var_param)]), max(df_results[paste0('true_', var_param)]))
  
  # LABELS
  if (FIXED_PARAM$r0){
    ylab <- expression(paste('Posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX = bquote(bold(paste(.(model) ~ "Model - ", italic(R[0]), " Inference")))
  } else {
    ylab <- paste0('Posterior mean of ', fixed_param)
    titleX = bquote(bold(paste(.(model) ~ "Model - " ~ .(fixed_param) ~ "Inference")))
  }
  
  #x label
  if (VAR_PARAM$r0) { #Var param should be x axis 
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
  } else {
    xlab <- paste0('True ', var_param)
  }
  
  # Create the plot for each subset
  for (i in 1:length(filter_list)) {
    
    df_subset <- filter_list[[i]]
    print(paste0('n row:', nrow(df_subset)))
    var_subset <- unlist(df_subset[paste0('true_', var_param)])
    mean_subset <- unlist(df_subset[paste0('mean_', fixed_param)])
    upper_ci_subset <- unlist(df_subset[paste0('upper_ci_', fixed_param)])
    lower_ci_subset <- unlist(df_subset[paste0('lower_ci_', fixed_param)])
    colour <- selected_colors[i]
    
    if (i == 1) {
      
      plot(var_subset, mean_subset, type = "p",
           main = "",
           xlab = xlab, ylab = ylab,
           xlim = x_lim,
           ylim = y_lim,
           col = colour, pch = 16, 
           #text.font = 4.0,
           cex = cex)
      title(main = list(titleX, cex = 1.2, font = 2.0))
      
    } else {
      
      points(var_subset, mean_subset, type = "p",
             xlab = xlab, ylab = ylab,
             ylim = c(0, max(var_subset, upper_ci_subset)),
             col = colour, pch = 16, cex = cex)
    }
    
    #Error bars 
    segments(var_subset, upper_ci_subset, var_subset,
             lower_ci_subset, lwd = 0.5, col = colour)
    # Points
    points(var_subset, mean_subset, type = "p",
           col = colour, pch = 16)
  }
  # TRUE line
  lines(var_total, fixed_total, col = 'black', lwd = 3)
  
  # Legend
  legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         inset=c(-0.32,0),
         col = c('black', selected_colors),
         lwd = rep(1.8, num_conds),
         lty = rep(1, num_conds), #c(1, 1),
         pch = c(NA, 19, 19, 19, 19, 19),
         text.font = 1.0,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}
