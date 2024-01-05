#PLOT INFERENCE RESULTS
'PLOT COMPUTED INFERENCE RESULTS FOR ALL MODELS '

PLOT_INFERENCE_RESULTS <- function(df_results, COMP_FOLDER, fig_num = '1',
                              num_days = 50, cex = 1.75, #1.25 
                              PDF = TRUE, GT = FALSE, GT_VAL = 20, inset = 0.46, #0.46 for non r0 
                              INCLUDE_INFS_5 = FALSE,
                              PRIORS = list(EXP = TRUE,
                                            GAMMA = FALSE, UNIF = FALSE, BETA = FALSE, GAMMA_B = FALSE),
                              FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE,
                                                alpha = FALSE, beta = FALSE),
                              FLAG_FILTER = list(end_day = FALSE,
                                                 tot_infs = TRUE),
                              FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                                                SSEB = FALSE, SSIB = FALSE)){
  
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
  par(mar=c(5.2, 4.8, 3.0, 19.0), xpd=TRUE) #Margins; bottom, left, top, right
  
  
  #DATA SUBSETS
  list_subset_data = SUBSET_DFS(df_results, filter_param, true_param, GT = GT, GT_VAL = GT_VAL)
  subset_df_list = list_subset_data$subset_df_list; legend_list = list_subset_data$legend_list
  selected_colors = list_subset_data$selected_colors; num_conds = list_subset_data$num_conds 
  
  #PLOT
  y_lim = c(0, max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
  x_lim = c(min(df_results[paste0('true_', true_param)]), max(df_results[paste0('true_', true_param)]))
  #y_lim = c(0,7) #c(0, max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
  
  #PRIORS
  prior_title = GET_PRIOR_TITLE(PRIORS)
  
  # LABELS
  if (FLAG_PARAM$r0) {
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
    ylab <- expression(paste('Estimated posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX =  bquote(paste(italic(R[0]) ~ " - " ~ .(model)))
    #titleX = bquote(bold(paste(italic(R[0])))) #, " Inference"))) #. Epidemic length: " ~ .(num_days) ~ ' days')))
    #titleX =  bquote(bold(paste(.(model) ~ "Model Inference - ", italic(R[0]))))
    
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
           cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = cex) #cex = cex text.font = 4.0,
      title(main = list(titleX, cex = 2.5, font = 3.0))
      
    } else {
      
      points(true_subset, mean_subset, type = "p",
             xlab = xlab, ylab = ylab,
             ylim = c(0, max(true_subset, upper_ci_subset)),
             col = colour, pch = 16,
             cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = cex)
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

#SUBSET THE DATASET
SUBSET_DFS <- function(df_results, filter_param, param,
                       GT=FALSE, GT_VAL = 20, num_conds = 4,
                       INCLUDE_INFS_5 = FALSE, FIXED = FALSE){
  
  'SUBSET DATAFRAME OF RESULTS'
  #Setup
  viridis_palette <- viridis(10)
  
  #SUBSET DATAFRAME OF RESULTS
  subset_df_list <- list(
    df10 = df_results[(df_results[[filter_param]] > 10) & (df_results[[filter_param]] <= 100), ],
    df100 = df_results[(df_results[[filter_param]] > 100) & (df_results[[filter_param]] <= 1000), ],
    df1k = df_results[(df_results[[filter_param]] > 1000) & (df_results[[filter_param]] <= 10000), ],
    df10k = df_results[(df_results[[filter_param]] > 10000) & (df_results[[filter_param]] <= 100000), ]) 
  
  #Legend list
  key = ' estimated mean, infections: ' #est infs
  if (INCLUDE_INFS_5){
    
    #ADD TO 
    subset_df_list$df5 = df_results[(df_results[[filter_param]] > 4) & (df_results[[filter_param]] <= 10), ]
    
    legend_list <- c(paste0(param, " True "),
                     paste0(param, " estimated mean, infections: [5, 10]"),
                     paste0(param, key, "[10, 100]"), 
                     paste0(param, key, "[100, 1k]"),
                     paste0(param, key, "[1k, 10k]"),
                     paste0(param, key, "[10k, 100k]"))
    
    num_conds = num_conds + 1
    selected_colors <- viridis_palette[c(9, 10, 8, 5)]
    
  } else if(GT) {
    
    legend_list <- c(paste0(param, " True "),
                     paste0(param, " estimated mean, infections: [", GT_VAL, ",100]"), 
                     paste0(param, key, "[100, 1k]"),
                     paste0(param, key, "[1k, 10k]"),
                     paste0(param, key, "[10k, 100k]"))
    
    selected_colors <- viridis_palette[c(10, 8, 5)] 
  } else {
    
    legend_list <- c(paste0(param, " True "),
                     paste0(param, " estimated mean, infections: [10, 100]"), 
                     paste0(param, key, "[100, 1k]"),
                     paste0(param, key, "[1k, 10k]"),
                     paste0(param, key, "[10k, 100k]"))
    
    selected_colors <- viridis_palette[c(10, 8, 5)]
    
  }
  
  #ADD MAGMA COLOURS 
  mag_cols = magma(5)
  colors_magma <- mag_cols[c(3, 4)]
  selected_colors = c(selected_colors, colors_magma)
  
  #Check for 100k
  if (!is.null(df_results[df_results[[filter_param]] > 100000, ])) {
    subset_df_list$df100k <- df_results[df_results[[filter_param]] > 100000, ]
    legend_list = c(legend_list, paste0(param, " estimated mean, infections > 100k"))
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
SIM_PERFORMANCE_R0 <- function(df_results){
  
  #Bias, MAE
  df_results$MAE = abs(df_results$mean_r0 - df_results$true_r0)
  df_results$bias = df_results$mean_r0 - df_results$true_r0
  num_runs = length(df_results$true_r0)
  
  #Results
  print(paste0('mean bias: ', round(mean(df_results$bias), 3)))
  print(paste0('MAE: ', round(mean(df_results$MAE), 3)))
  print(paste0('mean sd: ', round(mean(df_results$sd_r0), 3)))
  print(paste0('coverage: ', sum(df_results$coverage_r0)))
  print(paste0('% coverage: ', sum(df_results$coverage_r0)/num_runs))
  
  #return(df_results)
  
}

SIM_PERFORMANCE <- function(df_results){
  
  #Bias, MAE
  df_results$MAE = abs(df_results$mean_r0 - df_results$true_r0)
  df_results$bias = df_results$mean_r0 - df_results$true_r0
  num_runs = length(df_results$true_r0)
  
  #Results
  print(paste0('mean bias: ', round(mean(df_results$bias), 3)))
  print(paste0('MAE: ', round(mean(df_results$MAE), 3)))
  print(paste0('mean sd: ', round(mean(df_results$sd), 3)))
  print(paste0('coverage: ', sum(df_results$coverage)))
  print(paste0('% coverage: ', sum(df_results$coverage)/num_runs))
  
  #return(df_results)
  
}


#*****************************
#* PRIORS
#*****************************
PLOT_PRIOR_DIST <- function(PRIORS = list(EXP = FALSE, GAMMA = FALSE, UNIF = FALSE,
                                          BETA = FALSE, GAMMA_B = TRUE),
                            x_min = 0, x_max = 30, cex = 1.6){
  
  #Setup
  par(mfrow = c(2,3))
  x = seq(from = x_min, to = x_max, length = 500)
  #prior = names(PRIORS)[which(unlist(PRIORS))]
  
  #PRIORS
  if(PRIORS$EXP){
    print('a')
    y = dexp(x, 1)
    prior_title =  paste0('Exponential(1) Prior')
    
  } else if (PRIORS$GAMMA){
    prior_title =  paste0('Gamma(1, 5) Prior')
    y = dgamma(x, shape = 1, scale = 5)
    
  } else if (PRIORS$UNIF){
    prior_title =  paste0('Uniform(0, 10) Prior')
    y = dunif(x, 0, 10)
  } else if (PRIORS$BETA){
    prior_title =  paste0('Beta(2, 2) Prior')
    y = dbeta(x, 2, 2)
  } else if (PRIORS$GAMMA_B){
    prior_title =  paste0('1 + Gamma(3, 3) Prior')
    y = dgamma(x-1, shape = 3, scale = 3)
  }
  
  plot(x, y, type = 'l', lwd = 2, #col = 'orange',
       main = prior_title, ylab = 'Density',
       cex.lab=cex, cex.axis=cex, cex.main= cex, cex.sub=cex)
}

GET_PRIOR_TITLE <-function(PRIORS){
  
  #PRIORS
  if(PRIORS$EXP){
    prior_title =  'Exponential(1) Prior used'
  } else if (PRIORS$GAMMA){
    prior_title =  'Gamma(1, 5) Prior used'
  } else if (PRIORS$UNIF){
    prior_title =  'Uniform(0, 10) Prior used'
  } else if (PRIORS$BETA){
    prior_title =  'Beta(2, 2) Prior used'
  } else if (PRIORS$GAMMA_B){
    prior_title =  '1 + Gamma(3, 3) Prior used' 
  }
}