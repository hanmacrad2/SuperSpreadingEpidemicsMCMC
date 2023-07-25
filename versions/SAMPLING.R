#SUBSAMPLE OF DATA

vector = data_baseline

calculate_magnitude <- function(x) {
  abs(x)
}

total_magnitude <- sum(sapply(vector, calculate_magnitude))
target_sum <- 0.8 * total_magnitude

sampled_vector <- c()
current_sum <- 0
running_prob_sum <- 0

for (value in vector) {
  magnitude <- calculate_magnitude(value)
  probability <- magnitude / total_magnitude
  
  if ((running_prob_sum + probability) >= 0.8) {
    break
  }
  
  sampled_vector <- c(sampled_vector, value)
  current_sum <- current_sum + magnitude
  running_prob_sum <- running_prob_sum + probability
}

p_value <- running_prob_sum
p_value

#VERSION 1
vector <- c(2, 0, 0, 0, 1, 0, 2, 2, 1, 1, 0, 0, 0, 2, 1, 2, 1, 2, 5, 1, 3, 5, 2, 4, 2, 3, 4, 1, 4, 7, 4, 5, 6, 5, 4, 9, 8, 10, 5, 11, 9, 12, 12, 15, 13, 16, 20, 14, 17, 19)

p = 0.7426471
calculate_magnitude <- function(x) {
  abs(x)
}

total_magnitude <- sum(sapply(vector, calculate_magnitude))
target_sum <- 0.8 * total_magnitude

sampled_vector <- list()
current_sum <- 0

for (value in vector) {
  magnitude <- calculate_magnitude(value)
  probability <- magnitude / total_magnitude
  
  if (probability >= p) {
    sampled_vector <- c(sampled_vector, value)
    current_sum <- current_sum + magnitude
    
    if (current_sum >= target_sum) {
      break
    }
  }
}

if (current_sum > target_sum) {
  sampled_vector <- sampled_vector[-length(sampled_vector)]
}

sampled_vector