#MODEL CRITICISM: POSTERIOR PREDICTIVE CHECKS
library(purrr)

#FUNCTIONS
rdata <- function() {
  rnorm(num_days, 0, sd = 1)
}

rtraj <- function() {
  rnorm(n = num_days, sd = 2)
}

zigzag <- function(xs) {
  sum(abs(diff(xs)))
}

#PARAMS
num_days <- 20
num_reps <- 100

#DATA
true_data <- rdata()
post_samps <- lapply(1:num_reps, \(n) rtraj()) #n represents the current iteration value

#Zig-zag
apply(m, 1, zigzag)

post_pred_samp_zz <- post_samps |> map(zigzag) |> unlist()


m <- matrix(0, nrow = num_reps, ncol = num_days)
for (i in 1:num_reps) {
  m[i, ] <- as.numeric(post_samps[[i]])
}

upper_bounds <- apply(m, 2, quantile, probs = 0.975) |> unlist()

mean_est <- apply(m, 2, mean) |> unlist()

lower_bounds <- apply(m, 2, quantile, probs = 0.025) |> unlist()

#PLOTS 
par(mfrow = c(1, 2))
hist(post_pred_samp_zz, xlim = c(0, max(max(post_pred_samp_zz), zigzag(true_data) * 2)))
abline(v = quantile(post_pred_samp_zz, probs = c(0.025, 0.975)), col = 'red')
abline(v = zigzag(true_data))

plot(1:num_days, upper_bounds, type = 'l', ylim = c(-5, 5))
lines(1:num_days, rep(qnorm(p=0.025, sd = 2), num_days), type = 'l')
lines(1:num_days, mean_est, type = 'l')
lines(1:num_days, rep(qnorm(p=0.975, sd = 2), num_days), type = 'l')
lines(1:num_days, lower_bounds, type = 'l')
points(1:num_days, true_data, type = 'l', col = 'red')