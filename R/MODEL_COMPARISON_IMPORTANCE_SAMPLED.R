#Model Comparison - Importance Sampled

#LIBRARIES
library(mvtnorm)

#Trial
alpha = c(0.8, 0.81, 0.805, 0.9)
beta = c(0.02, 0.025, 0.03, 0.5)
gamma = c(10, 10.2, 10.1, 11)
mcmc_output = matrix(c(alpha, beta, gamma), ncol = 3)
#Mean of mcmc
means = colMeans(mcmc_output)

#DEFENSE MIXTURE
t_dist = rmvt(100, sigma = cov(mcmc_output), df = 3) + means
out = 0.95*dmvt(t_dist - means, sigma = cov(mcmc_output), df = 3, log = FALSE)


out#LOAD MCMC
