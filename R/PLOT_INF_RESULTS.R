#PLOT INFERENCE RESULTS

#R0
PLOT_CI_R0_FILTER <- function(df_results, title, cex = 0.8, 
                     FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = FALSE),
                      list_colors = list(col_true = 'black',
                  col_ci = 'orange', col_mean_gt ='red', col_mean_lt = 'grey')){
  
  #par(mfrow=c(1,1))
  
  #Filter
  if(FLAG_FILTER$SUM){
    df_gt = df_results[df_results$tot_infs > 3, ]
    df_lt = df_results[df_results$tot_infs  <= 3 , ] 
    leg_filter = 'Mean k, tot infs < 3'
  } else if (FLAG_FILTER$FINAL_DAY){
    df_gt = df_results[df_results$end_day > 0, ]
    df_lt = df_results[df_results$end_day  <= 0 , ] 
    leg_filter = 'Mean k, final day > 0'
  }

  #Plot Means (gt & lt)
  plot(df_gt$true_r0, df_gt$mean_r0, type = "p",
       main = title, #0.1 - 0.9. True(blk) Mean (red)' ),
       xlab = 'true r0', ylab = 'mean r0',
       ylim = c(0, max(df_gt$true_r0, df_gt$upper_ci_r0)),
       col = list_colors$col_mean_gt, pch = 16, cex = cex)
  points(df_lt$true_r0, df_lt$mean_r0, type = "p",
         col = list_colors$col_mean_lt, pch = 16)
  
  # Add error bars using segments()
  segments(df_gt$true_r0, df_gt$lower_ci_r0, df_gt$true_r0,
           df_gt$upper_ci_r0, lwd = 1, col = list_colors$col_ci) #col_ci)
  segments(df_lt$true_r0, df_lt$lower_ci_r0, df_lt$true_r0,
           df_lt$upper_ci_r0, lwd = 1, col = list_colors$col_mean_lt)
  
  #ADD AGAIN
  points(df_gt$true_r0, df_gt$mean_r0, type = "p",
         col = list_colors$col_mean_gt, pch = 16)
  points(df_lt$true_r0, df_lt$mean_r0, type = "p",
         col = list_colors$col_mean_lt, pch = 16)
  
  #TRUE
  lines(df_results$true_r0, df_results$true_r0, col = list_colors$col_true, lwd = 4)
  
  #Legend
  legend("topleft", legend = c("True k", "Mean k mcmc", "k mcmc CIs", 
                               leg_filter),
         col = c(list_colors$col_true, list_colors$col_mean_gt,
                 list_colors$col_ci, list_colors$col_mean_lt), lwd = c(2, 2, 2, 2),
         lty = c(1,1,1,1),
         pch = c(NA, 19, NA, 19), text.font = 2, bty = "n")
}

#1. Tot infs > 3
title = 'Baseline Model - R0 Inference'
pdf_file = "fig_baseline_sum_gt3.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_R0_FILTER(df_results, title,
                        FLAG_FILTER = list(FINAL_DAY = FALSE, SUM = TRUE))
dev.off()

#2. Final day > 0
title = 'Baseline Model - R0 Inference'
pdf_file = "fig_baseline_end_gt0.pdf" 
pdf(paste0(COMP_FOLDER, pdf_file), width = 6, height = 5)
PLOT_CI_R0_FILTER(df_results, title,
                  FLAG_FILTER = list(FINAL_DAY = TRUE, SUM = FALSE))
dev.off()