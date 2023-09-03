#************************
#* SSE INF 

#RUN_INFERENCE_SSE
df_sse_inf = RUN_INFERENCE_SSE(DATA_FOLDER) #, level = 0.90)
df_results = df_sse_inf

#********************************
#PLOT CIs
#********************************

cap = 0.03
plot(df_results$true_k, df_results$mean_k, type = "o",
     main = 'k Inference; SSE Model',
     xlab = 'true k', ylab = 'mean k',
     pch = 16, ylim = c(0, max(df_results$upper_ci_k)))

# Add error bars using segments()
segments(df_results$true_k, df_results$lower_ci_k, df_results$true_k, df_results$upper_ci_k, lwd = 2, col = "blue")

# Add horizontal lines/caps at the top of the vertical lines
segments(df_results$true_k - cap, 
         df_results$upper_ci_k, df_results$true_k + cap,
         df_results$upper_ci_k, lwd = 1, col = "blue")

segments(df_results$true_k - cap, df_results$lower_ci_k,
         df_results$true_k + cap,
         df_results$lower_ci_k, lwd = 1, col = "blue")

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


#*****************
#* RO PLOT
#* **************
#1. Add true r0 column


