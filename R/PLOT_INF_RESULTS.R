#PLOT INFERENCE RESULTS

#ALL PARAMS: NEW COLOUR
PLOT_CI_PARAMS <- function(df_results, title, cex = 0.8, 
                            FLAG_PARAM = list(r0 = FALSE, k = FALSE,  alpha = FALSE, gamma = FALSE),
                            FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = FALSE),
                            list_colors = list(col_true = 'black',
                            col_ci = 'orange', col_mean_gt ='orange',
                            col_mean_lt = 'grey')){
  
  #PARAM
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  #Colours
  viridis_palette <- viridis(10, alpha = 0.7)
  selected_colors <- viridis_palette[c(5, 9)]
  list_colors$col_ci = selected_colors[1]; list_colors$col_mean_gt = selected_colors[1]; 
  list_colors$col_mean_lt = selected_colors[2]
  
  #Filter
  if(FLAG_FILTER$SUM){
    df_gt = df_results[df_results$tot_infs > 3, ]
    df_lt = df_results[df_results$tot_infs  <= 3 , ] 
    leg_filter = paste0(true_param, ' Mean MCMC, total infs < 3')
  } else if (FLAG_FILTER$FINAL_DAY){
    df_gt = df_results[df_results$end_day > 0, ]
    df_lt = df_results[df_results$end_day  <= 0 , ] 
    leg_filter = paste0(true_param, ' Mean MCMC, epidemic dies out')
  }
  
  #EXTRACT
  true = unlist(df_gt[paste0('true_', true_param)])
  #true = unlist(df_gt[paste0(true_param)])
  mean = unlist(df_gt[paste0('mean_', true_param)])
  upper_ci = unlist(df_gt[paste0('upper_ci_', true_param)])
  lower_ci = unlist(df_gt[paste0('lower_ci_', true_param)])
  true_lt = unlist(df_lt[paste0('true_', true_param)])
  #true_lt = unlist(df_lt[paste0(true_param)])
  mean_lt = unlist(df_lt[paste0('mean_', true_param)])
  upper_ci_lt = unlist(df_lt[paste0('upper_ci_', true_param)])
  lower_ci_lt = unlist(df_lt[paste0('lower_ci_', true_param)])
  true_total = unlist(df_results[paste0('true_', true_param)])
  legend_list = c(paste0(true_param, " True "), paste0(true_param, " Mean mcmc"), #paste0(true_param, " MCMC CI "), 
                  leg_filter)
  #LABELS
  if(FLAG_PARAM$r0) {
    xlab = paste0('true R0'); ylab = paste0('mean R0') 
  } else {
    xlab = paste0('true ', true_param); ylab = paste0('mean', true_param)
  }
 
  #PLOTs
  plot(true, true, type = "p",
       main = title, #0.1 - 0.9. True(blk) Mean (red)' ),
       xlab = xlab, ylab = ylab,
       ylim = c(0, max(true, upper_ci)),
       col = list_colors$col_mean_gt, pch = 16, cex = cex)
  points(true_lt, mean_lt, type = "p",
         col = list_colors$col_mean_lt, pch = 16)
  
  # Add error bars using segments()
  segments(true, upper_ci, true,
           lower_ci, lwd = 0.5, col = list_colors$col_ci) #col_ci)
  segments(true_lt, lower_ci_lt, true_lt,
           upper_ci_lt, lwd = 0.5, col = list_colors$col_mean_lt)
  
  #ADD AGAIN
  points(true, mean, type = "p",
         col = list_colors$col_mean_gt, pch = 16)
  points(true_lt, mean_lt, type = "p",
         col = list_colors$col_mean_lt, pch = 16)
  
  #TRUE
  lines(true_total, true_total, col = list_colors$col_true, lwd = 4)
  
  #Legend
  legend("topleft", legend_list,
         col = c(list_colors$col_true, list_colors$col_mean_gt, #list_colors$col_ci
                 list_colors$col_mean_lt), lwd = c(1.5, 1.5, 1.5), # 2),
         lty = c(1,1,1),
         pch = c(NA, 19, 19), text.font = 1.5, bty = "n")
}


#FIXED PARAM
PLOT_CI_PARAM_FIXED <- function(df_results, title, cex = 0.8, 
                                FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE, gamma = FALSE),
                                FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = FALSE),
                                list_colors = list(col_true = 'black',
                                                   col_ci = 'orange', col_mean_gt ='orange',
                                                   col_mean_lt = 'grey')){
  
  'For fixed param only difference is seq along'
  #PARAM
  true_param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  #Colours
  viridis_palette <- viridis(10)
  selected_colors <- viridis_palette[c(5, 9)]
  
  #Filter
  if(FLAG_FILTER$SUM){
    df_gt = df_results[df_results$tot_infs > 3, ]
    df_lt = df_results[df_results$tot_infs  <= 3 , ] 
    leg_filter = paste0(true_param, ' Mean MCMC, total infs < 3')
  } else if (FLAG_FILTER$FINAL_DAY){
    df_gt = df_results[df_results$end_day > 0, ]
    df_lt = df_results[df_results$end_day  <= 0 , ] 
    leg_filter = paste0(true_param, ' Mean MCMC, epidemic dies out')
  }
  
  #EXTRACT
  true = unlist(df_gt[paste0('true_', true_param)])
  mean = unlist(df_gt[paste0('mean_', true_param)])
  upper_ci = unlist(df_gt[paste0('upper_ci_', true_param)])
  lower_ci = unlist(df_gt[paste0('lower_ci_', true_param)])
  true_lt = unlist(df_lt[paste0('true_', true_param)])
  mean_lt = unlist(df_lt[paste0('mean_', true_param)])
  upper_ci_lt = unlist(df_lt[paste0('upper_ci_', true_param)])
  lower_ci_lt = unlist(df_lt[paste0('lower_ci_', true_param)])
  true_total = unlist(df_results[paste0('true_', true_param)])
  legend_list = c(paste0(true_param, " True "), paste0(true_param, " Mean mcmc"), #paste0(true_param, " MCMC CI "), 
                  leg_filter)
  #LABELS
  if(FLAG_PARAM$r0) {
    xlab = paste0('true R0'); ylab = paste0('mean R0') 
  } else {
    xlab = paste0('true ', true_param); ylab = paste0('mean', true_param)
  }
  
  #PLOTs
  plot(seq_along(true), mean, type = "p",
       main = title, #0.1 - 0.9. True(blk) Mean (red)' ),
       xlab = 'iter', ylab = ylab,
       ylim = c(0, max(true, upper_ci)),
       col = list_colors$col_mean_gt, pch = 16, cex = cex)
  points(seq_along(true_lt), mean_lt, type = "p",
         col = list_colors$col_mean_lt, pch = 16)
  
  # Add error bars using segments()
  segments(seq_along(true), upper_ci, seq_along(true),
           lower_ci, lwd = 0.5, col = list_colors$col_ci) #col_ci)
  segments(seq_along(true_lt), lower_ci_lt, seq_along(true_lt),
           upper_ci_lt, lwd = 0.5, col = list_colors$col_mean_lt)
  
  #ADD AGAIN
  points(seq_along(true), mean, type = "p",
         col = list_colors$col_mean_gt, pch = 16)
  points(seq_along(true_lt), mean_lt, type = "p",
         col = list_colors$col_mean_lt, pch = 16)
  
  #TRUE
  lines(seq_along(true_total), true_total, col = list_colors$col_true, lwd = 4)
  
  #Legend
  legend("topleft", legend_list,
         col = c(list_colors$col_true, list_colors$col_mean_gt, #list_colors$col_ci
                 list_colors$col_mean_lt), lwd = c(1.5, 1.5, 1.5), # 2),
         lty = c(1,1,1),
         pch = c(NA, 19, 19), text.font = 1.5, bty = "n")
}


#******************
#* 1. BASELINE 
#*****************

#1. Tot infs > 3
title = 'Baseline Model - R0 Inference'
pdf_file = "3fig_baseline_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_base, title, FLAG_PARAM = list(r0 = TRUE, k = FALSE, 
                                                 alpha = FALSE, gamma = FALSE),
                        FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

#2. Final day > 0
title = 'Baseline Model - R0 Inference'
pdf_file = "2fig_baseline_end_gt0.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_base, title,  FLAG_PARAM = list(r0 = TRUE, k = FALSE),
                  FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()

#******************
#* 2. SSE 
#*****************
title = 'SSE Model - k Inference'
pdf_file = "fig_sse_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_sse, title, FLAG_PARAM = list(r0 = FALSE, k = TRUE),
                  FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

#2. Final day > 0
title = 'SSE Model - k Inference'
pdf_file = "fig_sse_end_gt0.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_sse, title,  FLAG_PARAM = list(r0 = FALSE, k = TRUE),
                  FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()

#RENAME
df_sse$true_k <- df_sse$k
df_sse$k <- NULL 

#***********
#FIXED SSE: R0
#***********
title = 'SSE Model - R0 Inference. k:0-1'
pdf_file = "fig_sse1_r0_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_sse1, title, FLAG_PARAM = list(r0 = TRUE, k = FALSE),
               FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

title = 'SSE Model - R0 Inference. k:0-1'
pdf_file = "fig_sse1_r0_dies.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_sse1, title, FLAG_PARAM = list(r0 = TRUE, k = FALSE),
                    FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()

#************************
#* SSEB
#* **********************
#R0
title = 'SSEB Model - R0 Inference'
pdf_file = "fig_sseb_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_sseb, title, FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                          alpha = FALSE, gamma = TRUE),
               FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

#2. Final day > 0
title = 'SSEB Model' # - R0 Inference'
pdf_file = "fig_sseb_end_gt0.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_sseb, title,  FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                            alpha = FALSE, gamma = TRUE),    
               FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()

#***********
#OTHER FIXED PARAMS SSEB
#***********
title = 'SSEB Model'
pdf_file = "fig_sse1_r0_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_sseb, title, 
                    FLAG_PARAM = list(r0 = FALSE, k = FALSE,
                                      alpha = FALSE, gamma = TRUE),
                    FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

title = 'SSE Model - R0 Inference. k:0-1'
pdf_file = "fig_sse1_r0_dies.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_sseb, title, FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                                                      gamma = TRUE),
                    FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()

#gamma
title = 'SSEB Model'
pdf_file = "fig_sse1_r0_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_sseb, title, 
                    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE, gamma = TRUE),
                    FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

title = 'SSE Model - R0 Inference. k:0-1'
pdf_file = "fig_sse1_r0_dies.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_sseb, title, FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE, gamma = TRUE),
                    FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()


#***********
#* SSI
#***********
title = 'SSI Model - k Inference'
pdf_file = "fig_ssi1_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_ssi1, title, FLAG_PARAM = list(r0 = FALSE, k = TRUE),
               FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

#2. Final day > 0
title = 'SSI Model - k Inference'
pdf_file = "fig_ssi1_end_gt0.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAMS(df_ssi1, title,  FLAG_PARAM = list(r0 = FALSE, k = TRUE),
               FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()

#***********
#FIXED SSI: R0
#***********
title = 'SSI Model - R0 Inference'
pdf_file = "fig_ssi_r0_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_ssi, title, FLAG_PARAM = list(r0 = TRUE, k = FALSE,
                                                     alpha = FALSE, gamma = FALSE),
                    FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

title = 'SSI Model - R0 Inference' #. k:0-1'
pdf_file = "fig_ssi_r0_dies.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_PARAM_FIXED(df_ssi, title, FLAG_PARAM = list(r0 = TRUE, k = FALSE,
                                                     alpha = FALSE, gamma = FALSE),
                    FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()