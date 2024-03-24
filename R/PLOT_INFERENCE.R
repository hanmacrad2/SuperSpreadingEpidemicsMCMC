#**********************************
#*
#PLOT INFERENCE RESULTS

#**********************************
PLOT_INFERENCE_RESULTS <- function(df_results, COMP_FOLDER, fig_num = '1',
                                  cex = 1.7, PDF = TRUE, GT_VAL = 30, inset = 0.485, #cex = 1.75
                                  PLOT_PRIOR_CI = TRUE,
                                  FLAG_PARAM = GET_PARAM(r0 = TRUE), 
                                  FLAG_MODEL = GET_FLAGS_MODELS(BASELINE = TRUE), FIXED = FALSE){
  
  'Plot computed Inference results for all models '
  #PARAMS 
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  true_total = unlist(df_results[paste0('true_', true_param)])
  model = names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  #FLAG_PRIOR = GET_PRIOR(EXP = FALSE)
  filter_param = 'tot_infs'
  
  #PDF
  if(PDF){
    plot_folder = paste0(COMP_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_', true_param, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = 13.0, height = 8.0) #13, 8
  }
  par(mar=c(5.2, 4.8, 3.0, 19.45), xpd=TRUE) #Margins; bottom, left, top, right
  
  
  #DATA SUBSETS
  list_subset_data = SUBSET_DFS(df_results, filter_param, true_param, FLAG_PARAM, GT_VAL = GT_VAL)
  subset_df_list = list_subset_data$subset_df_list; legend_list = list_subset_data$legend_list
  num_conds = list_subset_data$num_conds
  selected_colors =  GET_SELECTED_COLORS()
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
           cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = cex) #cex = cex text.font = 4.0,
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
  }
  
  #TRUE line
  lines(true_total, true_total, col = 'black', lwd = 3)
  
  #PLOT PRIOR_CI
  dark_gray =  gray(0.5)
  if(PLOT_PRIOR_CI){
    PLOT_PRIOR_CI(FLAG_PARAM, dark_gray)
  }
  
  #Legend
  GET_LEGEND_II(legend_list, prior_title, selected_colors, dark_gray, num_conds)
  # legend_list = c(legend_list, '', prior_title)
  # pch_list = rep(19, num_conds + 1)
  # legend('bottomright', #x = "topleft", y = "topleft", #"center", legend_list,
  #        legend_list,
  #        cex = 1.1,
  #        inset=c(-inset,0),
  #        col = c('black', selected_colors, 'white', dark_gray),
  #        lwd = rep(1.8, num_conds+2),
  #        lty = rep(1, num_conds+2), #c(1, 1),
  #        pch = c(NA, pch_list, NA, NA),
  #        #text.font = 1.8, #1.45
  #        bty = "n")
  
  if(PDF){
    dev.off()
  }
}


#LEGEND FUNCTION
GET_LEGEND_II <- function(legend_list, prior_title,
                          selected_colors, dark_gray, num_conds, inset = 0.485,
                       legend_location = 'bottomright'){
  
  #COLOURS ADD WHITE
  selected_colors = c(selected_colors[1], 'white',
                      selected_colors[2], 'white',
                      selected_colors[3], 'white',
                      selected_colors[4], 'white',
                      selected_colors[5], 'white',
                      dark_gray)
  
  #Legend
  legend_list = c(legend_list, '', prior_title)
  pch_list = rep(19, num_conds + 1)
  legend(legend_location, 
         legend_list,
         cex = 1.1,
         inset=c(-inset,0),
         col = c('black', selected_colors, 'white', dark_gray),
         lwd = rep(1.8, num_conds+2),
         lty = rep(1, num_conds+2), 
         pch = c(NA, pch_list, NA, NA),
         bty = "n")
}

#PLOT_PRIOR_CI
PLOT_PRIOR_CI <- function(FLAG_PARAM, dark_gray){
  
  #PLOT 95% CIS OF PRIORS
  x0 <- 0 # Position it at the start of the x-axis
  
  if(FLAG_PARAM$r0){
    y_start <- 0.025
    y_end <- 3.688
  } else if (FLAG_PARAM$k) {
    y_start <- 0.051
    y_end <- 0.693
  } else if (FLAG_PARAM$alpha | FLAG_PARAM$a){
    y_start <- 0.035
    y_end <- 0.965
  } else if (FLAG_PARAM$beta | FLAG_PARAM$b){
    y_start <- GET_GAMMA_CI()[1]
    y_end <- GET_GAMMA_CI()[2]
  }

  #PLOT VERTICAL LINE
  segments(x0, y_start, x0, y_end, col = dark_gray, lwd = 2)
  
}

#GET SELECTED COLOURS
GET_SELECTED_COLORS <-function(){
  
  #PALETTE
  viridis_palette <- viridis(10)
  selected_colors <- viridis_palette[c(10, 8, 5)] 
  
  #ADD MAGMA COLOURS 
  mag_cols = magma(5)
  colors_magma <- mag_cols[c(3, 4)]
  selected_colors = c(selected_colors, colors_magma)
  
  return(selected_colors)
}

#SUBSET THE DATASET
SUBSET_DFS <- function(df_results, filter_param, param, FLAG_PARAM,
                       GT = TRUE, GT_VAL = 30, num_conds = 5, #5
                       FIXED = FALSE) { #FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE, alpha = FALSE, beta = FALSE)){
  
  'SUBSET DATAFRAME OF RESULTS'
  #Setup
  #viridis_palette <- viridis(10)
  
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
  
  #selected_colors <- viridis_palette[c(10, 8, 5)] 
  
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
  

  #Check for 100k
  subset_df_list$df100k <- df_results[df_results[[filter_param]] > 100000, ]

  
  return(list(subset_df_list = subset_df_list, legend_list = legend_list, 
              num_conds = num_conds))
             # selected_colors = selected_colors, num_conds = num_conds))
}

#POSTERIOR + PRIOR
SCALE_PARAM <- function(vec_param){
  
  vec_scaled = vec_param/max(vec_param)
  
  return(vec_scaled)
}

#****************
#* POSTERIOR & PRIOR PLOTS + HISTOGRAMS
PLOT_HIST_PRIOR <- function(df_results, FLAG_PARAM, FLAGS_MODELS, MODEL_COLOR,
                            RESULTS_FOLDER, xlimits, ylimits, sim_val, 
                            n_repeats = 1000, main_font = 1.5, legend_location = 'topright',
                            inset = 0.2, cex = 1.8, bar_height = 20, true_height = 0.5, alpha = 0.2, 
                            border_alpha = NA, PDF = TRUE){
  
  #MODEL
  num_reps =  nrow(df_results); print(paste0('N reps: ', num_reps))
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))] 
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  COLOR_ALPHA = GET_COLOR_ALPHA(MODEL_COLOUR, alpha)
  COLOR_BODER = GET_COLOR_ALPHA(MODEL_COLOUR, 0.3)
  mcmc_vec1 = unlist(df_results[paste0(param, '_mcmc')][1,])
  
  #*********************
  
  if(PDF){
    
    #FOLDER
    plot_folder = paste0(RESULTS_FOLDER, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_', param, '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = 12.5, height = 7.0) #13, 8
    par(mar=c(5.2, 4.8, 3.0, 19.45), xpd=TRUE) #Margins; bottom, left, top, right
  }
 
  
  #HISTOGRAM
  #title = bquote(.(list_labels$main_inf), ~ '. N repeats = ' ~.(n_repeats))
  hist(mcmc_vec1, freq = FALSE,
       col = COLOR_ALPHA,
       xlim = xlimits, ylim = ylimits,
       xlab = list_labels$lab, 
       ylab = 'Estimated Posterior Density',
       border = NA,
       cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = 2.5,
       main = list(list_labels$main_inf, cex = 1.9, font = main_font)) 
  
  #HISTOGRAMS ADDITIONAL 

  for (i in c(2:n_repeats)){
    mcmc_vec = unlist(df_results[paste0(param, '_mcmc')][i,])
    hist(mcmc_vec, freq = FALSE,
         col = COLOR_ALPHA,
         xlim = xlimits, ylim = ylimits,
         xlab = list_labels$lab, 
         ylab = 'Estimated Posterior Density',
         border = border_alpha, #NA,
         cex.lab=cex, cex.axis=cex-0.3, cex.sub=cex-0.3, cex = 2.5,
         add = TRUE)
    
  }
  
  #PLOT PRIOR
  PLOT_PRIOR_DIST(FLAG_PARAM, xlimits = xlimits, ylimits = ylimits, alpha = 0.3)
  
  #PLOT TRUE
  max_y = ylimits[2]
  segments(sim_val, 0, sim_val, bar_height, col = 'black', lwd = 1.85) #2 #max_y-true_height
  
  #LEGEND
  prior_title = GET_PRIOR_TITLE(FLAG_PARAM)
  legend_list = c(list_labels$legend_posterior, prior_title,
                  paste0('True Simulated Value: ', sim_val))
  #legend_list = c('alpha', prior_title,
  #                paste0('True Simulated Value: ', sim_val))
  
  if(PDF){
    GET_LEGEND(legend_list, COLOR_ALPHA, legend_location, inset)
  }

  
  if(PDF){
    dev.off()
  }
  
}

#LEGEND FUNCTION
GET_LEGEND <- function(legend_list, COLOR_ALPHA,
                       legend_location = 'topright', inset = 0.25){
  
  #COLOUR
  COLOR_ALPHA = GET_COLOR_ALPHA(COLOR_ALPHA, alpha = 0.8)
  
  #Legend
  num_conds = length(legend_list)
  pch_list = rep(19, num_conds)
  legend(legend_location, #x = "topleft", y = "topleft", #"center", legend_list,
         legend_list,
         cex = 1.1,
         inset=c(-inset,0),
         col = c(COLOR_ALPHA, 'gray', 'black'),
         lwd = rep(3, num_conds-1), #c(rep(3, num_conds-1), 2),
         lty = rep(1, num_conds), #c(1, 1),
         #pch = pch_list, #c(NA, pch_list, NA, NA, NA),
         text.font = 1.8, #1.45
         bty = "n")
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

#SUBSET THE DATASET
GET_SUBSET_DATAFRAME <- function(df_results, filter_param, param, FLAG_PARAM,
                       GT = TRUE, GT_VAL = 30, num_conds = 5, #5
                       FIXED = FALSE) { #FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE, alpha = FALSE, beta = FALSE)){
  
  'SUBSET DATAFRAME OF RESULTS'
  
  #SUBSET DATAFRAME OF RESULTS
  subset_df_list <- list(
    df10 = df_results[(df_results[[filter_param]] > 10) & (df_results[[filter_param]] <= 100), ],
    df100 = df_results[(df_results[[filter_param]] > 100) & (df_results[[filter_param]] <= 1000), ],
    df1k = df_results[(df_results[[filter_param]] > 1000) & (df_results[[filter_param]] <= 10000), ],
    df10k = df_results[(df_results[[filter_param]] > 10000) & (df_results[[filter_param]] <= 100000), ]) 
  
  return(subset_df_list)
  
}
#   #Legend list
#   key = ' estimated mean, infections:' #est infs
#   
#   #LEGEND LIST - DEFAULT 
#   legend_list <- c(paste0(param, " True "),
#                    paste0(param, " estimated mean & 95% CI, "), 
#                    paste0("infections: [", GT_VAL, ",100]"), 
#                    #paste0(param, " estimated mean, infections: [", GT_VAL, ",100]"), 
#                    #paste0(param, key, "[100, 1k]"),
#                    paste0(param, key, "[1k, 10k]"),
#                    paste0(param, key, "[10k, 100k]"),
#                    paste0(param, " estimated mean, infections > 100k")) 
#   
#   selected_colors <- viridis_palette[c(10, 8, 5)] 
#   
#   #R0 ALPHA OR BETA: EXPRESSIONS 
#   if(FLAG_PARAM$r0) {
#     legend_list <- c(expression(paste(R[0], " True ")),
#                      bquote(R[0] ~ .(key) ~ "[" * .(GT_VAL) * ",100]"),
#                      #bquote(R[0]~.(key)~"["~.(GT_VAL)~",100]"),
#                      bquote(R[0]~.(key)~"[100,1k]"),
#                      bquote(R[0]~.(key)~"[1k,10k]"),
#                      bquote(R[0]~.(key)~"[10k,100k]"),
#                      bquote(R[0]~"estimated mean, infections > 100k")) 
#     
#   } else if(FLAG_PARAM$alpha){
#     legend_list <- c(expression(paste(alpha, " True ")),
#                      bquote(alpha~.(key)~"["~.(GT_VAL) ~ ",100]"),
#                      bquote(alpha~.(key)~"[100, 1k]"),
#                      bquote(alpha~.(key)~"[1k, 10k]"),
#                      bquote(alpha~.(key)~"[10k, 100k]"),
#                      bquote(alpha~"estimated mean, infections > 100k")) 
#     
#   } else if (FLAG_PARAM$beta){
#     legend_list <- c(expression(paste(beta, " True ")),
#                      bquote(beta~.(key)~"["~.(GT_VAL)~",100]"),
#                      bquote(beta~.(key)~"[100,1k]"),
#                      bquote(beta~.(key)~"[1k,10k]"),
#                      bquote(beta~.(key)~"[10k,100k]"),
#                      bquote(beta~"estimated mean, infections > 100k")) 
#   }
#   
#   #ADD MAGMA COLOURS 
#   mag_cols = magma(5)
#   colors_magma <- mag_cols[c(3, 4)]
#   selected_colors = c(selected_colors, colors_magma)
#   
#   #Check for 100k
#   subset_df_list$df100k <- df_results[df_results[[filter_param]] > 100000, ]
#   
#   
#   return(list(subset_df_list = subset_df_list, legend_list = legend_list, 
#               selected_colors = selected_colors, num_conds = num_conds))
# }

#GET LEGEND LISTS
GET_LEGEND_LIST <- function(param, FLAG_PARAM,
                       GT = TRUE, GT_VAL = 30) { #FLAG_PARAM = list(r0 = FALSE, k = FALSE, kappa = FALSE, alpha = FALSE, beta = FALSE)){

  #LEGEND LIST
  lineI =  paste0(param, " estimated mean & 95% CIs, ")
  legend_list <- c(paste0(param, " True "),
                   lineI, 
                   paste0("infections: [", GT_VAL, ",100]"), 
                   lineI,
                   "infections: [100, 1k]",
                   lineI,
                   "infections: [1k, 10k]",
                   lineI,
                   "infections: [10k, 100k]",
                   lineI,
                   paste0("infections > 100k")) 

  
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
  
}