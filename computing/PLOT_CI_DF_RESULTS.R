#PLOT CI
library(gplots)
library(ggplot2)

#**********************
#FOLDERS
FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/computing/mcmc/"

#BASELINE
folder = '~/PhD_Warwick/Project_Epidemic_Modelling/Project_details/Computing/mcmc/baseline/'
df_file = 'base_mcmc_2023-08-10_17-50-05.rds'
df_results = readRDS(paste0(folder, df_file))

#SSE
FOLDER = paste0(FOLDER, 'sse/')
df_file = 'sse_mcmc_2023-08-11_18-45-13.rds'
df_results = readRDS(paste0(FOLDER, df_file))

#PLOT
df_results['true_k'] = rep(0.1, length(df_results$true_r0))
df_results$mean_r0_ci <- (df_results$lower_ci_r0 + df_results$upper_ci_r0) / 2
df_results$mean_k_ci <- (df_results$lower_ci_k + df_results$upper_ci_k) / 2

#PLOT (MAKE FUNCTION)
pchX = 16; lwdX = 1
#plotCI(seq_along(df_results[['true_r0']]), df_results[['mean_ci']],

#R0
plotCI(x = df_results[['true_r0']], y = df_results[['true_r0']],
       ui = df_results$upper_ci_r0, li = df_results$lower_ci_r0,
       xlab = 'r0 true', ylab = 'r0 true',
       gap = 0.0, #KEY! :D
       main = 'R0 Inference, SSE model. 101 simulations (unseen).',
      lwd = 1, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))

points(df_results[['true_r0']], df_results[['true_r0']], col = 'red', lwd = lwdX, pch = 16)
lines(df_results[['true_r0']], col = 'red', lwd = lwdX)

#k
plotCI(x = df_results[['lower_ci_k']], y = df_results[['lower_ci_k']],
       ui = df_results$upper_ci_k, li = df_results$lower_ci_k,
       xlab = 'k true', ylab = 'k true',
       gap = 0.0, #KEY! :D
       main = 'k Inference, SSE model. 101 simulations (unseen).')
#       lwd = 1, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))

abline(h = 0.1, col = 'red')
points(df_results[['true_k']], df_results[['true_r0']], col = 'red', lwd = lwdX, pch = 16)
lines(df_results[['true_k']], col = 'red', lwd = lwdX)

#k II
plotCI(x = seq_along(df_results[['lower_ci_k']]), y = df_results[['mean_k_ci']],
       ui = df_results$upper_ci_k, li = df_results$lower_ci_k,
       xlab = 'k true', ylab = 'k true',
       gap = 0.0, #KEY! :D
       main = 'k Inference, SSE model. 101 simulations (unseen).')
#       lwd = 1, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))

abline(h = 0.1, col = 'red')
points(df_results[['true_k']], df_results[['true_r0']], col = 'red', lwd = lwdX, pch = 16)
lines(df_results[['true_k']], col = 'red', lwd = lwdX)

#GGPLOT
p <- ggplot(df_results, aes(x = true_)) +
  geom_line(aes(y = true_r0), color = "red") +
  geom_errorbar(aes(ymin = lower_ci_r0, ymax = upper_ci_r0), width = 0.1) +
  labs(x = "True r0", y = "Value") +
  theme_minimal()

print(p)

