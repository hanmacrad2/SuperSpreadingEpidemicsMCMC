#APPLICATION OF MODEL EVIDENCE ESTIMATE (VIA IMPORTANCE SAMPLING)
library(compositions)

#1. LOAD IN MCMC 
mcmc = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/rjmcmc', i, '.rds' ))
PLOT_SSEB_RJMCMC(data_sseb1, rj_sse10, n_mcmc)

#2. CALCULATE PHATs








#MULTIVARIATE LOG NORMAL
#NOT MULTI-VARIATE
#multivariable log normal
MyVar <- matrix(c(
  0.2,0.1,0.0,
  0.1,0.2,0.0,
  0.0,0.0,0.2),byrow=TRUE,nrow=3)
MyMean <- c(1,1,2)

#Samples
rlnorm.rplus(100,log(MyMean),MyVar)

plot.ts(rlnorm.rplus(100,log(MyMean),MyVar))

plot(rlnorm.rplus(100,log(MyMean),MyVar))
plot(rnorm.aplus(100,MyMean,MyVar))
x <- rnorm.aplus(5,MyMean,MyVar)
dnorm.aplus(x,MyMean,MyVar)
dlnorm.rplus(x,log(MyMean),MyVar)