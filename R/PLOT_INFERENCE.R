#PLOT INFERENCE RESULTS
'PLOT COMPUTED INFERENCE RESULTS FOR ALL MODELS '

PLOT_INFERENCE_RESULTS <- function(df_results, COMP_FOLDER,fig_num = '1',
                              num_days = 50, cex = 1.25, 
                              PDF = TRUE, GT_20 = TRUE, INCLUDE_INFS_5 = FALSE,
                              PRIORS = list(EXP = TRUE,
                                            GAMMA = FALSE, UNIF = FALSE),
                              FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                                alpha = FALSE, gamma = FALSE),
                              FLAG_FILTER = list(end_day = FALSE,
                                                 tot_infs = TRUE),
                              FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                                SSEB = FALSE, SSIB = FALSE)){
  
  print(GT_20)
  #PARAMS 
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  true_total = unlist(df_results[paste0('true_', true_param)])
  filter_param = names(FLAG_FILTER)[which(unlist(FLAG_FILTER))]
  model =  names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  
  #PLOT
  plot_folder = paste0(COMP_FOLDER, '/plots/')
  create_folder(plot_folder)
  
  if(PDF){
    pdf_file = paste0(model, '_', true_param, '_', fig_num, '.pdf') #'Fig_', 
    pdf(paste0(plot_folder, pdf_file), width = 13.0, height = 8.0)
  }
  par(mar=c(4.9, 4.6, 3.0, 17.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  
  #DATA SUBSETS
  list_subset_data = SUBSET_DFS(df_results, filter_param, true_param, GT_20 = GT_20)
  subset_df_list = list_subset_data$subset_df_list; legend_list = list_subset_data$legend_list
  selected_colors = list_subset_data$selected_colors; num_conds = list_subset_data$num_conds 
  
  #PLOT
  #y_lim = c(0, max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
  y_lim = c(0, 2+max(true_total, max(df_results[paste0('mean_', true_param)]))) 
  x_lim = c(min(df_results[paste0('true_', true_param)]), max(df_results[paste0('true_', true_param)]))
  
  #PRIORS
  prior_title = GET_PRIOR_TITLE(PRIORS)
  
  # LABELS
  if (FLAG_PARAM$r0) {
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
    ylab <- expression(paste('Estimated posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX = bquote(bold(paste(.(model) ~ "Model Inference - ", italic(R[0])))) #, " Inference"))) #. Epidemic length: " ~ .(num_days) ~ ' days')))
    #titleX = bquote(paste(titleX, prior_title))
    
  } else {
    xlab <- paste0('True ', true_param)
    ylab <- paste0('Estimated posterior mean of ', true_param)
    titleX = bquote(bold(paste(.(model) ~ "Model Inference - " ~ .(true_param) ))) #~
    #"Inference"))) #. Epidemic length: " ~ .(num_days) ~ ' days')))
  }
  
  # Create the plot for each subset
  for (i in 1:length(subset_df_list)) {
    
    df_subset <- subset_df_list[[i]]
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
           col = colour, pch = 16, 
           cex.lab=cex, cex.axis=cex, cex.sub=cex) #cex = cex text.font = 4.0,
      title(main = list(titleX, cex = 1.8, font = 2.0))
      
    } else {
      
      points(true_subset, mean_subset, type = "p",
             xlab = xlab, ylab = ylab,
             ylim = c(0, max(true_subset, upper_ci_subset)),
             col = colour, pch = 16, cex = cex)
    }
    
    #ERROR BARS using segments()
    segments(true_subset, upper_ci_subset, true_subset,
             lower_ci_subset, lwd = 0.5, col = colour)
    
    #POINTS
    points(true_subset, mean_subset, type = "p",
           col = colour, pch = 16)
  }
  
  #TRUE line
  lines(true_total, true_total, col = 'black', lwd = 3)
  
  #Legend
  legend_list = c(legend_list, prior_title)
  pch_list = rep(19, num_conds)
  legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         inset=c(-0.40,0),
         col = c('black', selected_colors, 'grey'),
         lwd = rep(1.8, num_conds+1),
         lty = rep(1, num_conds+1), #c(1, 1),
         pch = c(NA, pch_list, NA),
         text.font = 1.2,
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}

#SUBSET THE DATASET
SUBSET_DFS <- function(df_results, filter_param, true_param, GT_20 = GT_20, num_conds = 4,
                       INCLUDE_INFS_5 = FALSE ){
  
  'SUBSET DATAFRAME OF RESULTS'
  print(GT_20)
  #Setup
  viridis_palette <- viridis(10)
  
  #SUBSET DATAFRAME OF RESULTS
  subset_df_list <- list(
    df10 = df_results[(df_results[[filter_param]] > 10) & (df_results[[filter_param]] <= 100), ],
    df100 = df_results[(df_results[[filter_param]] > 100) & (df_results[[filter_param]] <= 1000), ],
    df1k = df_results[(df_results[[filter_param]] > 1000) & (df_results[[filter_param]] <= 10000), ],
    df10k = df_results[(df_results[[filter_param]] > 10000) & (df_results[[filter_param]] <= 100000), ]) 
  
  #Legend list
  key = ' est mean, infs: '
  if (INCLUDE_INFS_5){
    
    #ADD TO 
    subset_df_list$df5 = df_results[(df_results[[filter_param]] > 4) & (df_results[[filter_param]] <= 10), ]
    
    legend_list <- c(paste0(true_param, " True "),
                     paste0(true_param, " estimated mean, infs: [5, 10]"),
                     paste0(true_param, key, "[10, 100]"), 
                     paste0(true_param, key, "[100, 1k]"),
                     paste0(true_param, key, "[1k, 10k]"),
                     paste0(true_param, key, "[10k, 100k]"))
    
    num_conds = num_conds + 1
    selected_colors <- viridis_palette[c(9, 10, 8, 5)]
    
  } else if(GT_20) {
    
    legend_list <- c(paste0(true_param, " True "),
                     paste0(true_param, " estimated mean, infections: [20,100]"), 
                     paste0(true_param, key, "[100, 1k]"),
                     paste0(true_param, key, "[1k, 10k]"),
                     paste0(true_param, key, "[10k, 100k]"))
    
    selected_colors <- viridis_palette[c(10, 8, 5)] 
  } else {
    
    legend_list <- c(paste0(true_param, " True "),
                     paste0(true_param, " estimated mean, infections: [10, 100]"), 
                     paste0(true_param, key, "[100, 1k]"),
                     paste0(true_param, key, "[1k, 10k]"),
                     paste0(true_param, key, "[10k, 100k]"))
    
    selected_colors <- viridis_palette[c(10, 8, 5)]
    
  }
  
  #ADD MAGMA COLOURS 
  mag_cols = magma(5)
  colors_magma <- mag_cols[c(3, 4)]
  selected_colors = c(selected_colors, colors_magma)
  
  #Check for 100k
  if (!is.null(df_results[df_results[[filter_param]] > 100000, ])) {
    subset_df_list$df100k <- df_results[df_results[[filter_param]] > 100000, ]
    legend_list = c(legend_list, paste0(true_param, " est mean, infs > 100k"))
    num_conds = num_conds + 1
  }
  
  return(list(subset_df_list = subset_df_list, legend_list = legend_list, 
              selected_colors = selected_colors, num_conds = num_conds))
}

#***********************
#*
#* PLOT PERFORMANCE & INFERENCE RESULTS
#* 
#* **********************
SIM_PERFORMANCE <- function(df_results){
  
  #Bias, MAE
  df_results$MAE = abs(df_results$mean_r0 - df_results$true_r0)
  df_results$bias = df_results$mean_r0 - df_results$true_r0
  
  #Results
  print(paste0('mean bias: ', round(mean(df_results$bias), 3)))
  print(paste0('MAE: ', round(mean(df_results$MAE), 3)))
  print(paste0('mean sd: ', round(mean(df_results$sd), 3)))
  print(paste0('coverage: ', sum(df_results$coverage)))
  print(paste0('% coverage: ', sum(df_results$coverage)/1000))
  
  #return(df_results)
  
}