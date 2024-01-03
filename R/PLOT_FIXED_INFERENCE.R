#********************************************************
#* PLOT FIXED PARAMETERS OF INFERENCE COMPUTATION RUNS
#********************************************************

PLOT_FIXED_INFERENCE <- function(df_results, COMP_FOLDER, fig_num = '1',
                                   num_days = 50, cex = 1.25, 
                                   PDF = TRUE, GT = FALSE, GT_VAL = 20, inset = 0.43, #0.46 for non r0 
                                   INCLUDE_INFS_5 = FALSE,
                                   PRIORS = list(EXP = TRUE,
                                                 GAMMA = FALSE, UNIF = FALSE, BETA = FALSE, GAMMA_B = FALSE),
                                   FIXED_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE,
                                                     alpha = FALSE, beta = FALSE, b = FALSE),
                                   VAR_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE,
                                                   alpha = FALSE, beta = FALSE, b = FALSE),
                                   FLAG_FILTER = list(end_day = FALSE,
                                                      tot_infs = TRUE),
                                   FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                                     SSEB = FALSE, SSIB = FALSE)){
  
  #PARAMS 
  var_param =  names(VAR_PARAM)[which(unlist(VAR_PARAM))]
  fixed_param = names(FIXED_PARAM)[which(unlist(FIXED_PARAM))]
  fixed_total = unlist(df_results[paste0('true_', fixed_param)])
  var_total = unlist(df_results[paste0('true_', var_param)])
  filter_param = names(FLAG_FILTER)[which(unlist(FLAG_FILTER))]
  model =  names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  
  #PLOT
  plot_folder = paste0(COMP_FOLDER, '/plots/')
  create_folder(plot_folder)
  
  if(PDF){
    pdf_file = paste0(model, '_var_', var_param, '_fixed_', fixed_param, '_', fig_num, '.pdf') #'Fig_', 
    pdf(paste0(plot_folder, pdf_file), width = 13.0, height = 8.0)
  }
  par(mar=c(4.9, 4.6, 3.0, 19.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  
  #DATA SUBSETS
  list_subset_data = SUBSET_DFS(df_results, filter_param, fixed_param, GT = GT, GT_VAL = GT_VAL)
  subset_df_list = list_subset_data$subset_df_list; legend_list = list_subset_data$legend_list
  selected_colors = list_subset_data$selected_colors; num_conds = list_subset_data$num_conds 
  
  #PLOT
  y_lim = c(0, max(fixed_total, max(df_results[paste0('upper_ci_', fixed_param)]))) 
  x_lim = c(min(df_results[paste0('true_', var_param)]), max(df_results[paste0('true_', var_param)]))
  
  #PRIORS
  prior_title = GET_PRIOR_TITLE(PRIORS)
  
  # LABELS
  if (FIXED_PARAM$r0){
    ylab <- expression(paste('Estimated posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX = bquote(bold(paste(.(model) ~ "Model Inference - ", italic(R[0]))))
  } else {
    ylab <- paste0('Estimated posterior mean of ', fixed_param)
    titleX = bquote(bold(paste(.(model) ~ "Model Inference - " ~ .(fixed_param) )))
  }
  
  #x label
  if (VAR_PARAM$r0) { #Var param should be x axis 
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
  } else {
    xlab <- paste0('True ', var_param)
  }
  
  # Create the plot for each subset
  for (i in 1:length(subset_df_list)) {
    
    df_subset <- subset_df_list[[i]]
    print(paste0('n row:', nrow(df_subset)))
    var_subset <- unlist(df_subset[paste0('true_', var_param)])
    mean_subset <- unlist(df_subset[paste0('mean_', fixed_param)])
    upper_ci_subset <- unlist(df_subset[paste0('upper_ci_', fixed_param)])
    lower_ci_subset <- unlist(df_subset[paste0('lower_ci_', fixed_param)])
    colour <- selected_colors[i]
    
    if (i == 1) {
      
      plot(var_subset, mean_subset, type = "p",
           main = "", xlab = xlab, ylab = ylab,
           xlim = x_lim, ylim = y_lim,
           col = colour, pch = 16, 
           cex.lab=cex, cex.axis=cex, cex.sub=cex) #cex = cex text.font = 4.0,
      title(main = list(titleX, cex = 1.8, font = 2.0))
      
    } else {
      
      points(var_subset, mean_subset, type = "p",
             xlab = xlab, ylab = ylab,
             ylim = c(0, max(var_subset, upper_ci_subset)),
             col = colour, pch = 16, cex = cex)
    }
    
    #ERROR BARS using segments()
    segments(var_subset, upper_ci_subset, var_subset,
             lower_ci_subset, lwd = 0.5, col = colour)
    
    #POINTS
    points(var_subset, mean_subset, type = "p",
           col = colour, pch = 16)
  }
  
  #TRUE line
  lines(var_total, fixed_total, col = 'black', lwd = 3)
  
  #Legend
  legend_list = c(legend_list, '', prior_title)
  pch_list = rep(19, num_conds)
  legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         inset=c(-inset,0),
         col = c('black', selected_colors, 'white', 'white'),
         lwd = rep(1.8, num_conds+2),
         lty = rep(1, num_conds+2), #c(1, 1),
         pch = c(NA, pch_list, NA, NA),
         text.font = 1.2,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}