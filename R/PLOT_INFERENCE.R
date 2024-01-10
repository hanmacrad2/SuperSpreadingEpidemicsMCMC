#PLOT INFERENCE RESULTS
'PLOT COMPUTED INFERENCE RESULTS FOR ALL MODELS '

PLOT_INFERENCE_RESULTS <- function(df_results, COMP_FOLDER, fig_num = '1',
                              num_days = 50, cex = 1.75, #1.25 
                              PDF = TRUE, GT_VAL = 20, inset = 0.486, #0.46 for non r0 
                              PRIORS = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE, BETA = FALSE, GAMMA_B = FALSE),
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
    pdf(paste0(plot_folder, pdf_file), width = 13.0, height = 8.0) #13, 8
  }
  par(mar=c(5.2, 4.8, 3.0, 19.45), xpd=TRUE) #Margins; bottom, left, top, right
  
  
  #DATA SUBSETS
  list_subset_data = SUBSET_DFS(df_results, filter_param, true_param, GT_VAL = GT_VAL, FLAG_PARAM = FLAG_PARAM)
  subset_df_list = list_subset_data$subset_df_list; legend_list = list_subset_data$legend_list
  selected_colors = list_subset_data$selected_colors; num_conds = list_subset_data$num_conds
  
  #PRIORS
  prior_title = GET_PRIOR_TITLE(PRIORS)
  
  # LABELS
  ylab = 'Estimated posterior mean of '
  if (FLAG_PARAM$r0) {
    xlab <- expression(paste('True R'[0])) #paste0(expression(paste('true R'[0])))  #
    ylab <- expression(paste('Estimated posterior mean of R'[0])) #paste0(expression(paste('mean R'[0]))) 
    titleX =  bquote(paste(italic(R[0]) ~ " - " ~ .(model)))
    #titleX = bquote(bold(paste(italic(R[0])))) #, " Inference"))) #. Epidemic length: " ~ .(num_days) ~ ' days')))
    #titleX =  bquote(bold(paste(.(model) ~ "Model Inference - ", italic(R[0]))))
    
  } else if (FLAG_PARAM$alpha){
    titleX =  bquote(paste(italic(alpha) ~ " - " ~ .(model)))
    xlab = bquote(paste("True " ~ italic(alpha)))
    ylab = bquote(paste(.(ylab) ~ italic(alpha)))
    
  } else if (FLAG_PARAM$beta) {
    titleX =  bquote(paste(italic(beta) ~ " - " ~ .(model)))
    xlab = bquote(paste("True " ~ italic(beta)))
    ylab = bquote(paste(.(ylab) ~ italic(beta)))
    
  } else {
    xlab <- paste0('True ', true_param)
    ylab <- paste0('Estimated posterior mean of ', true_param)
    titleX =  bquote(paste(.(true_param) ~ " - " ~ .(model))) #.(true_param)
    #titleX = bquote(bold(paste(.(model) ~ "Model Inference - " ~ .(true_param) ))) #~
    #"Inference"))) #. Epidemic length: " ~ .(num_days) ~ ' days')))
  }
  
  #LIMITS
  x_lim = c(min(df_results[paste0('true_', true_param)]), max(df_results[paste0('true_', true_param)]))
  
  if(FLAG_PARAM$k){
    y_lim = c(0, 2.0) #1.0; if 0.01-> 0.2s
  } else {
    y_lim = c(0, max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
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
           cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = 2.5) #cex = cex text.font = 4.0,
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
         cex = 1.1,
         inset=c(-inset,0),
         col = c('black', selected_colors, 'white', 'white'),
         lwd = rep(1.8, num_conds+2),
         lty = rep(1, num_conds+2), #c(1, 1),
         pch = c(NA, pch_list, NA, NA),
         #text.font = 1.8, #1.45
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}

#SUBSET THE DATASET
SUBSET_DFS <- function(df_results, filter_param, param,
                       GT=TRUE, GT_VAL = 10, num_conds = 5, 
                       FIXED = FALSE, FLAG_PARAM = FLAG_PARAM) { #FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE, alpha = FALSE, beta = FALSE)){
  
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
  key = 'estimated mean, infections:' #est infs
  
  #LEGEND LIST - DEFAULT 
  legend_list <- c(paste0(param, " True "),
                   paste0(param, " estimated mean, infections: [", GT_VAL, ",100]"), 
                   paste0(param, key, "[100, 1k]"),
                   paste0(param, key, "[1k, 10k]"),
                   paste0(param, key, "[10k, 100k]"))
  
  selected_colors <- viridis_palette[c(10, 8, 5)] 
  
  #R0 ALPHA OR BETA: EXPRESSIONS 
  if(FLAG_PARAM$r0) {
    legend_list <- c(expression(paste(R[0], " True ")),
                     bquote(R[0] ~ .(key) ~ "[" * .(GT_VAL) * ",100]"),
                     #bquote(R[0]~.(key)~"["~.(GT_VAL)~",100]"),
                     bquote(R[0]~.(key)~"[100,1k]"),
                     bquote(R[0]~.(key)~"[1k,10k]"),
                     bquote(R[0]~.(key)~"[10k,100k]"),
                     bquote(R[0]~"estimated mean, infections > 100k")) 
    
  } else if(FLAG_PARAM$alpha){
    legend_list <- c(expression(paste(alpha, " True ")),
                     bquote(alpha~.(key)~"["~.(GT_VAL) ~ ",100]"),
                     bquote(alpha~.(key)~"[100, 1k]"),
                     bquote(alpha~.(key)~"[1k, 10k]"),
                     bquote(alpha~.(key)~"[10k, 100k]"),
                     bquote(alpha~"estimated mean, infections > 100k")) 
    
  } else if (FLAG_PARAM$beta){
    legend_list <- c(expression(paste(beta, " True ")),
                     bquote(beta~.(key)~"["~.(GT_VAL)~",100]"),
                     bquote(beta~.(key)~"[100,1k]"),
                     bquote(beta~.(key)~"[1k,10k]"),
                     bquote(beta~.(key)~"[10k,100k]"),
                     bquote(beta~"estimated mean, infections > 100k")) 
  }

  #ADD MAGMA COLOURS 
  mag_cols = magma(5)
  colors_magma <- mag_cols[c(3, 4)]
  selected_colors = c(selected_colors, colors_magma)
  
  #Check for 100k
  subset_df_list$df100k <- df_results[df_results[[filter_param]] > 100000, ]
  # if (!is.null(df_results[df_results[[filter_param]] > 100000, ])) {
  #   subset_df_list$df100k <- df_results[df_results[[filter_param]] > 100000, ]
  #   legend_list = c(legend_list, paste0(param, " estimated mean, infections > 100k"))
  #   num_conds = num_conds + 1
  # }
  
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
#* PRIORS                       #FLAG_PARAM = list(r0 = TRUE, k = FALSE)
#*****************************
PLOT_PRIOR_DIST <- function(PRIORS = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE,
                                          BETA = FALSE, GAMMA_B = FALSE), cex = 2.65){
  
  #Setup
  par(mfrow = c(2,3))
  
  par(mar=c(5.2, 5.2, 5.2, 5.2), xpd=TRUE)
  #prior = names(PRIORS)[which(unlist(PRIORS))]
  
  #PRIORS
  if(PRIORS$EXP){
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    y = dexp(x, 1)
    x_label = expression(R[0])
    prior_title =  bquote("Prior; " ~ p[R[0]] ~ "= Exponential(1)")
    #prior_title =  paste0('Exponential(1) Prior')
    
  } else if (PRIORS$GAMMA){
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    prior_title =  paste0('Gamma(1, 5) Prior')
    y = dgamma(x, shape = 1, scale = 5)
    
  } else if (PRIORS$UNIF){
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    prior_title =  paste0('Uniform(0, 10) Prior')
    y = dunif(x, 0, 10)
  } else if (PRIORS$BETA){
    
    x_min = 0; x_max = 1
    x = seq(from = x_min, to = x_max, length = 500)
    #prior_title =  paste0('Beta(2, 2) Prior')
    #x_label = expression(beta)
    prior_title =  bquote("Prior; " ~ p[alpha] ~ "= Beta(2, 2)")
    y = dbeta(x, 2, 2)
    
  } else if (PRIORS$GAMMA_B){
    x_min = 0; x_max = 30
    x = seq(from = x_min, to = x_max, length = 500)
    x_label = expression(beta)
    prior_title =  bquote("Prior; " ~ p[beta] ~ "= 1 + Gamma(3, 3)")
    y = dgamma(x-1, shape = 3, scale = 3)
    #prior_title =  paste0('1 + Gamma(3, 3) Prior')
  }
  
  plot(x, y, type = 'l', lwd = 2, #col = 'orange',
       main = prior_title, 
       #xlab = x_label, 
       ylab = 'Density',
       cex.lab=cex, cex.axis=cex-0.3, cex.main= cex, cex.sub=cex-0.3)
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