#INSPECT PRIOR PLOTS
par(mfrow = c(1,1))
#EXP
# Grid of X-axis values
x <- seq(0, 10, 0.1)

# lambda = 0.1
rate_1 = 1
rate_2 = 0.1
rate_3 = 10
plot(x, dexp(x, rate = rate_1), type = "l", ylab = "", lwd = 2, col = "blue",
     main = 'Exponential priors', xaxt = "n", xlim = c(0,10))
# lambda = 1
lines(x, dexp(x, rate = rate_2), col = "red", lty = 1, lwd = 2)
lines(x, dexp(x, rate = rate_3), col = "green", lty = 1, lwd = 2)

# Adding a legend
legend("topright", c(expression(paste(, lambda)), "0.1", "1"),
       lty = c(0, 1, 1), col = c("blue", "red"), box.lty = 0, lwd = 2)

#GAMMA
shape <- 2; scale = 1
mean <- 1.6 #1.6
scale <- shape * mean

# Generate x-values for the plot
x <- seq(0, 10, length = 100)

# Calculate the PDF values using dgamma()
pdf_values <- dgamma(x-1, shape, scale)

pdf_values <- dbeta(x/1.2, 1, 2)/1.2
# Plot the gamma distribution
plot(x, pdf_values, type = "l", xlab = "x", ylab = "Density", main = "Gamma Distribution with Shape and Scale Parameters")
