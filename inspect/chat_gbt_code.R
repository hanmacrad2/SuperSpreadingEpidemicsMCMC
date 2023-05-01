#CHAT GBT CODE

my_matrix <- matrix(c(1, -2, 3, -4, 5, -6, 7, -8, 9), nrow = 3)

# Use apply function to count negative values in each column
neg_count <- apply(my_matrix, 2, function(x) sum(x < 0))

# Print the number of negative values in each column
neg_count
