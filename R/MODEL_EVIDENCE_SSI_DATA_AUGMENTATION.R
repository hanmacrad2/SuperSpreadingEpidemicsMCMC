#Model Evidence for SSI model with data augmentation

#Density function (mass function) for the multi-nomial dirichlet model

MULTI_DIRICHLET <- function(data_vector, shrink_factor = 0.1){
  
  data_vector = shrink_factor*data_vector
  term1 = factorial(sum(data_vector)-1)
  term2 = prod(factorial(data_vector-1))
  term3 = prod(data_vector^(data_vector-1))
  
  prob = (term1/term2)*term3
  
  return(prob)
}

#TOY EXAMPLE
data_ssi = c(200, 600, 200)

MULTI_DIRICHLET(data_ssi)

MULTI_DIRICHLET <- function(data_vector, shrink_factor = 0.1) {
  data_vector <- shrink_factor * data_vector
  
  log_term1 <- lgamma(sum(data_vector))
  log_term2 <- sum(lgamma(data_vector))
  log_term3 <- sum((data_vector - 1) * log(data_vector))
  
  log_prob <- log_term1 - log_term2 + log_term3
  prob <- exp(log_prob)
  
  return(prob)
}

# Example usage
data_ssi <- c(200, 600, 200)
prob <- MULTI_DIRICHLET(data_ssi)
print(prob)
