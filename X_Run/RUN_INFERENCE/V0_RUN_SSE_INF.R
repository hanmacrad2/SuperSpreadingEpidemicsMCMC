#************************
#* IMPLEMENT
#***********************
library(gplots)
library(ggplot2)

#PARAMS
lwd = 1

#FOLDER
folder_type = 'data_sse'
CURRENT_FOLDER = GET_FOLDER_TIME_STAMP(folder_type, array_index)
create_folder(CURRENT_FOLDER)

#RUN_INFERENCE_SSE
df_sse_inf = RUN_INFERENCE_SSE(DATA_FOLDER) #, level = 0.90)
df_results = df_sse_inf


#Dataframe columns
df_sse_inf$true_k <- df_sse_inf$true_r0
df_sse_inf$true_r0 <- NULL
df_sse_inf <- df_sse_inf[, c("true_k", "mean_k", "lower_ci_k", 'upper_ci_k',
                             'mean_r0',  "lower_ci_r0", "upper_ci_r0")]
df_results = df_sse_inf

#PLOT k CI
plotCI(x = df_results[['true_k']], y = df_results[['true_k']],
       ui = df_results$upper_ci_k, li = df_results$lower_ci_k,
       xlab = 'k true', ylab = 'k true',
       gap = 0.0, #KEY! :D
       main = 'k Inference, SSE model. 10 simulations.',
       lwd = 1, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))

points(df_results[['true_k']], df_results[['true_k']], col = 'red', lwd = lwd, pch = 16)
lines(df_results[['true_k']], df_results[['true_k']], col = 'red', lwd = lwd)


#GGPLOT2
p <- ggplot(df_results, aes(x = true_k)) +
  geom_line(aes(y = true_k), color = "red", size = 1.5) +
  geom_errorbar(aes(ymin = lower_ci_k, ymax = upper_ci_k), width = 0.1, size = 0.9) +  # Adjust the size value
  labs(x = "True k", y = "k") +
  ggtitle("k inference, SSE model  Mean, 95% CIs. (10 simulations)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

print(p)

#Plot
plot(x = df_results[['true_k']], y = df_results[['mean_k']], type = 'p')

#********************************
#PLOT OPTION 3
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


#*****************
#* RO PLOT
#* **************
#1. Add true r0 column


