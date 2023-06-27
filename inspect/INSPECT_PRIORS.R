#INSPECT PRIOR PLOTS

shape <- 2
mean <- 1.6 #1.6
scale <- shape * mean

# Generate x-values for the plot
x <- seq(0, 10, length = 100)

# Calculate the PDF values using dgamma()
pdf_values <- dgamma(x, shape, scale)

# Plot the gamma distribution
plot(x, pdf_values, type = "l", xlab = "x", ylab = "Density", main = "Gamma Distribution with Shape and Scale Parameters")