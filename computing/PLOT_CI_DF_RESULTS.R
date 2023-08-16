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
#df_file = 'sse_mcmc_2023-08-11_18-45-13.rds'
df_file = 'sse_mcmc_2023-08-14_17-13-02.rds'
df_file = 'sse_mcmc_2023-08-15_16-24-44.rds'
df_file = 'sse_mcmc_2023-08-16_11-46-37.rds'
df_results = readRDS(paste0(FOLDER, df_file))

#PLOT
df_results['true_k'] = rep(0.1, length(df_results$true_r0))
#df_results$mean_r0_ci <- (df_results$lower_ci_r0 + df_results$upper_ci_r0) / 2
#df_results$mean_k_ci <- (df_results$lower_ci_k + df_results$upper_ci_k) / 2

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
plotCI(x = seq_along(df_results[['lower_ci_k']]), y = df_results[['mean_k']],
       ui = df_results$upper_ci_k, li = df_results$lower_ci_k,
       xlab = 'iter', ylab = 'k',
       gap = 0.0, #KEY! :D
       main = 'k Inference, SSE model. k true: 0.1. Mean, 95% CIs.') #101 simulations (unseen).')
#       lwd = 1, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))

abline(h = 0.6, col = 'red', lwd = 3)
abline(h = 0.4, col = 'red', lwd = 3)
abline(h = 0.1, col = 'red', lwd = 3)
points(df_results[['true_k']], df_results[['true_r0']], col = 'red', lwd = lwdX, pch = 16)
lines(df_results[['true_k']], col = 'red', lwd = lwdX, type = 'o')

#PLOT
par(mfrow = c(1,1))
p <- ggplot(df_results, aes(x = true_r0)) +
  geom_line(aes(y = true_r0), color = "red", size = 1.5) +
  geom_errorbar(aes(ymin = lower_ci_r0, ymax = upper_ci_r0), width = 0.1, size = 0.9) +  # Adjust the size value
  labs(x = "True r0", y = "R0") +
  ggtitle("R0 inference, SSE model  Mean, 95% CIs. (101 unseen simulations)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

print(p)
plot(p)
