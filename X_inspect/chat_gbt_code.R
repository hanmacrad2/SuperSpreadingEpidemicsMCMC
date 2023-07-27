#CHAT GBT CODE

#NEGATIVE COUNT IN A COLUMN
my_matrix <- matrix(c(1, -2, Inf, -4, 5, -6, -Inf, -8, 0), nrow = 3)

# Use apply function to count negative values in each column
neg_count <- apply(my_matrix, 2, function(x) sum(x < 0))

zero_count <- apply(my_matrix, 2, function(x) sum(x == 0))

inf_count <- apply(my_matrix, 2, function(x) sum(is.infinite(x)))
inf_count



#PLOT WITH SUBTITLE
# create some sample data
x <- rnorm(100)

# create a histogram of the data with a main title and a subtitle
hist(x, main = "Histogram of x", sub = "Sample data", cex.sub = 0.8)
