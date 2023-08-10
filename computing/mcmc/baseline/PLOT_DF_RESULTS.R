#PLOT CI
library(gplots)
library(ggplot2)

#FOLDERS
folder_baseline = '~/PhD_Warwick/Project_Epidemic_Modelling/Project_details/Computing/mcmc/baseline/'
df_file = 'base_mcmc_2023-08-10_17-50-05.rds'
df_results = readRDS(paste0(folder_baseline, df_file))

#PLOT
df_results$mean_ci <- (df_results$lower_ci_r0 + df_results$upper_ci_r0) / 2

#PLOT
pchX = 16; lwdX = 1
#plotCI(seq_along(df_results[['true_r0']]), df_results[['mean_ci']],
plotCI(df_results[['true_r0']], df_results[['mean_ci']],
       ui = df_results$upper_ci_r0, li = df_results$lower_ci_r0,
       xlab = 'r0 true', ylab = 'r0 true',
       gap = 0.0, #KEY! :D
       main = 'R0 Inference, Baseline model',
      lwd = 1, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))

points(df_results[['true_r0']], df_results[['true_r0']], col = 'red', lwd = lwdX, pch = 16)
lines(df_results[['true_r0']], col = 'red', lwd = lwdX)

#GGPLOT
p <- ggplot(df_results, aes(x = true_r0)) +
  geom_line(aes(y = true_r0), color = "red") +
  geom_errorbar(aes(ymin = lower_ci_r0, ymax = upper_ci_r0), width = 0.1) +
  labs(x = "True r0", y = "Value") +
  theme_minimal()

print(p)