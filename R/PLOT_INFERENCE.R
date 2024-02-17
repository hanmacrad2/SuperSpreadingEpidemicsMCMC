#**********************************
#*
#PLOT INFERENCE RESULTS

#**********************************
PLOT_INFERENCE_RESULTS <- function(df_results, COMP_FOLDER, fig_num = '1',
                                  cex = 1.75, PDF = TRUE, GT_VAL = 10, inset = 0.485, 
                                   FLAG_PARAM = GET_PARAM(r0 = TRUE), 
                                  FLAG_MODEL = GET_FLAGS_MODELS(BASELINE = TRUE)){
  
  'Plot computed Inference results for all models '
  #PARAMS 
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  true_total = unlist(df_results[paste0('true_', true_param)])
  model = names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  #FLAG_PRIOR = GET_PRIOR(EXP = FALSE)
  filter_param = 'tot_infs'
  
  #PLOT
  plot_folder = paste0(COMP_FOLDER, '/plots/')
  create_folder(plot_folder)
  
  #PDF
  if(PDF){
    pdf_file = paste0(model, '_', true_param, '_', fig_num, '.pdf') #'Fig_', 
    pdf(paste0(plot_folder, pdf_file), width = 13.0, height = 8.0) #13, 8
  }
  par(mar=c(5.2, 4.8, 3.0, 19.45), xpd=TRUE) #Margins; bottom, left, top, right
  
  
  #DATA SUBSETS
  list_subset_data = SUBSET_DFS(df_results, filter_param, true_param, FLAG_PARAM, GT_VAL = GT_VAL)
  subset_df_list = list_subset_data$subset_df_list; legend_list = list_subset_data$legend_list
  selected_colors = list_subset_data$selected_colors; num_conds = list_subset_data$num_conds
  df_results = subset(df_results, tot_infs > GT_VAL)
  
  #PRIORS
  prior_title = GET_PRIOR_TITLE(FLAG_PARAM)
  
  #LIMITS
  x_lim = c(min(df_results[paste0('true_', true_param)]), max(df_results[paste0('true_', true_param)]))
  y_lim = c(min(true_total, min(df_results[paste0('lower_ci_', true_param)])),
            max(true_total, max(df_results[paste0('upper_ci_', true_param)]))) 
  
  #PLOT EACH SUBSET (DATA TOTAL)
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
           main = "", 
           xlab = list_labels$xlab, ylab = list_labels$ylab,
           xlim = x_lim, ylim = y_lim,
           col = colour, pch = 16, 
           cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = 2.5) #cex = cex text.font = 4.0,
      title(main = list(list_labels$main_inf, cex = 2.5, font = 3.0))
      
    } else {
      
      points(true_subset, mean_subset, type = "p",
             xlab = list_labels$xlab, ylab = list_labels$ylab,
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
         pch = c(NA, pch_list, NA, NA, NA),
         #text.font = 1.8, #1.45
         bty = "n")
  
  if(PDF){
    dev.off()
  }
}

#SUBSET THE DATASET
SUBSET_DFS <- function(df_results, filter_param, param, FLAG_PARAM,
                       GT = TRUE, GT_VAL = 10, num_conds = 5, #5
                       FIXED = FALSE) { #FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE, alpha = FALSE, beta = FALSE)){
  
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
  key = ' estimated mean, infections:' #est infs
  
  #LEGEND LIST - DEFAULT 
  legend_list <- c(paste0(param, " True "),
                   paste0(param, " estimated mean, infections: [", GT_VAL, ",100]"), 
                   paste0(param, key, "[100, 1k]"),
                   paste0(param, key, "[1k, 10k]"),
                   paste0(param, key, "[10k, 100k]"),
                   paste0(param, " estimated mean, infections > 100k")) 
  
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

  
  return(list(subset_df_list = subset_df_list, legend_list = legend_list, 
              selected_colors = selected_colors, num_conds = num_conds))
}

#POSTERIOR + PRIOR
SCALE_PARAM <- function(vec_param){
  
  vec_scaled = vec_param/max(vec_param)
  
  return(vec_scaled)
}

#****************
#* POSTERIOR & PRIOR PLOTS 
PLOT_POSTERIOR_PRIOR <- function(df_results, FLAG_PARAM, FLAGS_MODELS, MODEL_COLOR,
                                 xlimits, RESULTS_FOLDER, sim_val, main_font = 1.5,
                                 i = 2, cex = 1.8, alpha = 0.2, PDF = TRUE){
  
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))] 
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  COLOR_ALPHA = GET_COLOR_ALPHA(MODEL_COLOUR, alpha)
  limits = list(xlim = xlimits, ylim = c(0,1)) #GET_LIMITS(FLAG_PARAM, model)
  
  #PLOT
  plot_folder = paste0(RESULTS_FOLDER, '/plots/')
  create_folder(plot_folder)
  
  if(PDF){
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_', param, '_', time_stamp, '.pdf') #'Fig_', 
    pdf(paste0(plot_folder, pdf_file), width = 13.0, height = 8.0) #13, 8
  }
  par(mar=c(5.2, 4.8, 3.0, 19.45), xpd=TRUE) #Margins; bottom, left, top, right
  
  #1. DENSITY EXTRACT 
  density_mcmc = density(unlist(df_results[paste0(param, '_mcmc')][2,]))
  scaled_dx = density_mcmc$y/max(density_mcmc$y)
  #limits = list(xlim = c((min(density_mcmc$x)- 0.2), (max(density_mcmc$x) + 0.2)), ylim = c(0,1))
  
  plot(density_mcmc$x, scaled_dx,
       type = 'l', 
       col = COLOR_ALPHA,
       #col =  col.alpha(MODEL_COLOR, alpha = 0.15), 
       xlim = limits$xlim,
       xlab = list_labels$lab, 
       ylab = 'Estimated Posterior Density',
       cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = 2.5) 
  title(main = list(list_labels$main_inf, cex = 1.9, font = main_font))
  
  #POLYGON
  x = density_mcmc$x; y = scaled_dx
  polygon(c(x, rev(x)), c(y, rep(0, length(y))),
          col = COLOR_ALPHA, border = NA)
  
  #2. ADDITIONAL PLOTS
  for (i in c(3:50)){
    
    
    density_mcmc = density(unlist(df_results[paste0(param, '_mcmc')][i,]))
    scaled_dx = density_mcmc$y/max(density_mcmc$y)
    lines(density_mcmc$x, scaled_dx,
          type = 'l', col = MODEL_COLOR)
    #POLYGON FILL 
    x = density_mcmc$x; y = scaled_dx
    polygon(c(x, rev(x)), c(y, rep(0, length(y))),
            col = COLOR_ALPHA, border = NA)
  }
  
  #3. PRIOR
  #prior_xlim = c(limits$xlim[1] - 0.2, limits$xlim[2] + 0.2)
  prior_xlim = c(limits$xlim[1], limits$xlim[2])
  PLOT_PRIOR_DIST(FLAG_PARAM, xlimits = prior_xlim, alpha = 0.3)
  # x = seq(0, 6, by = 0.01)
  # dx1 = dexp(x, rate = 1)
  # dx2 = SCALE_PARAM(dx1)
  # lines(x, dx2, col = "gray", lty = "dashed", lwd = 2)
  
  #PLOT TRUE
  segments(sim_val, 0, sim_val, 1, col = 'black', lwd = 2)
  
  #LEGEND
  prior_title = GET_PRIOR_TITLE(FLAG_PARAM)
  legend_list = c(paste0('Estimated Posteriors of ', param), prior_title)
  GET_LEGEND(legend_list, COLOR_ALPHA)
  
  if(PDF){
    dev.off()
  }
  
}

#GET ALPHA COLOUR
GET_COLOR_ALPHA <- function(MODEL_COLOUR, alpha = 0.2){
  
  rgba_color <- col2rgb(MODEL_COLOR, alpha = TRUE)
  r = rgba_color[,]['red'][1]/255
  g = rgba_color[,]['green'][1]/255
  b = rgba_color[,]['blue'][1]/255
  
  COLOR_ALPHA = rgb(r, g, b, alpha = alpha)
  
  return(COLOR_ALPHA)
}

GET_LEGEND <- function(legend_list, COLOR_ALPHA){
  
  #COLOUR
  COLOR_ALPHA = GET_COLOR_ALPHA(COLOR_ALPHA, alpha = 0.8)
  
  #Legend
  num_conds = length(legend_list)
  pch_list = rep(19, num_conds)
  legend('topright', #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         cex = 1.1,
         #inset=c(-inset,0),
         col = c(COLOR_ALPHA, 'gray'),
         lwd = rep(3, num_conds),
         lty = rep(1, num_conds), #c(1, 1),
         #pch = pch_list, #c(NA, pch_list, NA, NA, NA),
         text.font = 1.8, #1.45
         bty = "n")
}