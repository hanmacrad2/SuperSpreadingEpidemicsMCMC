#PLOT INFERENCE RESULTS
PLOT_CI_PARAM_OTHER <- function(df_results, COMP_FOLDER, fig_num = '1',
                           cex = 1.0,
                           PDF = FALSE, 
                           FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                             alpha = FALSE, gamma = FALSE),
                           FLAG_FILTER = list(end_day = FALSE,
                                              tot_infs = FALSE,
                                              max_eta_x = TRUE), 
                           FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = TRUE,
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
  viridis_palette <- viridis(10) #inferno(10) #viridis, turbo
  selected_colors <- viridis_palette #[c(3, 4, 1, 8, 5, 7, 6, 10, 2, 9)] #[1:9] #[c(9, 10, 6, 8, 5)]
  num_conds = length(selected_colors)
  
  # Define different filters and their corresponding labels
  filter_list <- list(
    df1 = df_results[df_results[[filter_param]] < 0.5, ],
    df2 = df_results[(df_results[[filter_param]] >= 0.5) & (df_results[[filter_param]] < 1), ],
    df3 = df_results[(df_results[[filter_param]] >= 1) & (df_results[[filter_param]] < 2), ],
    df4 = df_results[(df_results[[filter_param]] >= 2) & (df_results[[filter_param]] <= 4), ],
    df5 = df_results[(df_results[[filter_param]] >= 5) & (df_results[[filter_param]] <= 10), ],
    df7 = df_results[(df_results[[filter_param]] > 10) & (df_results[[filter_param]] <= 15), ],
    df8 = df_results[(df_results[[filter_param]] > 15) & (df_results[[filter_param]] <= 20), ],
    df9 = df_results[(df_results[[filter_param]] > 20) & (df_results[[filter_param]] <= 25), ],
    df10 = df_results[df_results[[filter_param]] > 30, ]
  )
  #browser()
  
  # Create a list of legends for each subset
  #if(FLAG_FILTER$tot_infs){
  legend_list <- c(paste0(true_param, " True "),
                   paste0(true_param, " post. mean, max eta/x < 0.5"),
                   paste0(true_param, " post. mean, max eta/x: [0.5, 1]"),
                   paste0(true_param, " post. mean, max eta/x: [1, 2]"),
                   paste0(true_param, " post. mean, max eta/x: [2, 4]"),
                   paste0(true_param, " post. mean, max eta/x: [5, 10]"),
                   paste0(true_param, " post. mean, max eta/x: [10, 15]"),
                   paste0(true_param, " post. mean, max eta/x: [15, 20]"),
                   paste0(true_param, " post. mean, max eta/x: [20, 25]"),
                   paste0(true_param, " post. mean, max eta/x: > 30"))
                     
  #}
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
         pch = c(NA, rep(19, num_conds)), #c(NA, 19, 19, 19, 19, 19, 19, 10, ),
         text.font = 1.0,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}