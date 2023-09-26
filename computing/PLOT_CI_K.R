#PLOT CI k

#DATAFRAME
df_sse1 = df_results
file2 = 'sse_mcmc_2023-08-30_19-40-19.rds'
file2 = 'sse_mcmc_2023-08-31_17-02-05.rds'
df_sse2 = readRDS(paste0(COMP_FOLDER, file2))
df_results = df_sse2

#********************************
#PLOT CIs
#********************************
k_range = '1:10'
R0 = 1.5
cap = 0.009
max_infs = 300
col_ci = 'orange'; col_mean = 'red'
plot(df_results$k, df_results$mean_k, type = "p",
     main = paste0('SSE Model. k Inference. R0 = ', R0, ', k = ', k_range, '. Tot infs > ', max_infs), #0.1 - 0.9. True(blk) Mean (red)' ),
     xlab = 'true k', ylab = 'mean k',
     ylim = c(0, max(df_results$k, df_results$upper_ci_k)),
     col = col_mean, pch = 16)

# Add error bars using segments()
segments(df_results$k, df_results$lower_ci_k, df_results$k,
         df_results$upper_ci_k, lwd = 1, col = col_ci)

# Add horizontal lines/caps at the top of the vertical lines
segments(df_results$k - cap, 
         df_results$upper_ci_k, df_results$k + cap,
         df_results$upper_ci_k, lwd = 1, col = col_ci)

segments(df_results$k - cap, df_results$lower_ci_k,
         df_results$k + cap,
         df_results$lower_ci_k, lwd = 1, col = col_ci)

#ADD AGAIN
points(df_results$k, df_results$mean_k, type = "p",
       main = 'k Inference; SSE Model',
       xlab = 'true k', ylab = 'mean k',
       col = col_mean,
       pch = 16, ylim = c(0, max(df_results$upper_ci_k)))

#TRUE
col_true = 'black'
lines(df_results$k, df_results$k, col = col_true, lwd = 4)

#LEGEND
legend("topleft", legend = c("True k", "Mean k mcmc"),
       col = c(col_true, col_mean), lwd = c(3, 2), pch = c(NA, 19))

#********************************************************************************

#Axis ticks
ticks <- seq(0.1, 2, by = 0.1)
#labels <- as.character(ticks)

# Customize x-axis ticks
axis(1, at = ticks, labels = ticks)

points(df_results[['true_k']], df_results[['true_k']], col = 'red', lwd = lwd, pch = 16)
lines(df_results[['true_k']], df_results[['true_k']], col = 'red', lwd = lwd)

#axis(1, at = df_results$true_k, labels = df_results$true_k)
# Add legend
legend("topright", legend = c("True k", "Mean k"),
       col = c("red", "blue"), lwd = c(1, 2), pch = c(16, NA))


#*******************************

#********************************
#PLOT CIs
#********************************

cap = 0.03
plot(df_results$true_k, df_results$mean_k, type = "oo",
     main = 'k Inference; SSE Model, 90% CIs',
     xlab = 'true k', ylab = 'mean k',
     pch = 16, ylim = c(0, max(df_results$upper_ci_k_90)))

# Add error bars using segments()
segments(df_results$true_k, df_results$lower_ci_k_90, df_results$true_k, df_results$upper_ci_k_90, lwd = 2, col = "blue")

# Add horizontal lines/caps at the top of the vertical lines
segments(df_results$true_k - cap, 
         df_results$upper_ci_k_90, df_results$true_k + cap,
         df_results$upper_ci_k_90, lwd = 1, col = "blue")

segments(df_results$true_k - cap, df_results$lower_ci_k_90,
         df_results$true_k + cap,
         df_results$lower_ci_k_90, lwd = 1, col = "blue")

#Axis ticks
ticks <- seq(0.1, 2, by = 0.1)
#labels <- as.character(ticks)

# Customize x-axis ticks
axis(1, at = ticks, labels = ticks)

points(df_results[['true_k']], df_results[['true_k']], col = 'red', lwd = lwd, pch = 16)
lines(df_results[['true_k']], df_results[['true_k']], col = 'red', lwd = lwd)

#axis(1, at = df_results$true_k, labels = df_results$true_k)
# Add legend
legend("topright", legend = c("True k", "Mean k"),
       col = c("red", "blue"), lwd = c(1, 2), pch = c(16, NA))

legend("topright", legend = c("True k", "Mean k mcmc"),
       col = c('red', "blue"), lwd = c(1, 2), pch = c(16, NA))

#Load and save
file_name = 'df_sse_k_01_09_range_2023-08-30.rds'
df2 = readRDS(paste0(OUTER_FOLDER, file_name))
