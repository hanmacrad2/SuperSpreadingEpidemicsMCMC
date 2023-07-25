#INSPECT PRIOR PLOTS

#EXP
# Grid of X-axis values
x <- seq(0, 8, 0.1)

# lambda = 0.1
plot(x, dexp(x, rate = 1), type = "l", ylab = "", lwd = 2, col = "blue")
# lambda = 1
lines(x, dexp(x, rate = 0.1), col = "red", lty = 1, lwd = 2)

# Adding a legend
legend("topright", c(expression(paste(, lambda)), "2", "1"),
       lty = c(0, 1, 1), col = c("blue", "red"), box.lty = 0, lwd = 2)

#GAMMA
shape <- 2
mean <- 1.6 #1.6
scale <- shape * mean

# Generate x-values for the plot
x <- seq(0, 10, length = 100)

# Calculate the PDF values using dgamma()
pdf_values <- dgamma(x, shape, scale)

# Plot the gamma distribution
plot(x, pdf_values, type = "l", xlab = "x", ylab = "Density", main = "Gamma Distribution with Shape and Scale Parameters")