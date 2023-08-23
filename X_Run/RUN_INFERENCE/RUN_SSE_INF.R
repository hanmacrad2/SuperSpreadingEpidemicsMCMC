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
df_sse_inf = RUN_INFERENCE_SSE(DATA_FOLDER)

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


#*****************
#* RO PLOT
#* **************

