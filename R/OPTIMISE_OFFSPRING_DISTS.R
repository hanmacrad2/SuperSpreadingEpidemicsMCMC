#COST_FUNCTION
COST_FUNCTION <- function(list_params, sec_cases) {
  
  num_offspring = 1000
  x = 0:num_offspring
  z <- GET_OFFSPRING_SSIB(x, list_params)
  freq_z = round(z*length(x))
  
  sqrt(mean((sec_cases - freq_z)^2)) #RMSE 
  #sum((sec_cases - freq_z)^2) # Sum of squared differences
}

# Initial guess for the parameters
initial_guess <- c(0.5, 0.5, 10.0)

# Use optim to find the best parameters
fit_params <- optim(initial_guess, COST_FUNCTION, sec_cases = sec_cases)

result <- optim(initial_guess, COST_FUNCTION)

# Print the result
print(result)

#FREQUENCY
total_cases = sum(sec_cases_real)

# Convert probabilities to frequencies by scaling them up
freq <- (z*total_cases)


#OTHER
#Optim
poi_log_like <- function(r0_opt, y){
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda = r0_opt*sum(y[1:t-1]*rev(prob_infect[1:t-1]))
    logl = logl + y[t]*log(lambda) - lambda
    
  }
  
  return(-logl)
  
}

optim(1, poi_log_like, y = y)

result = optim(1, poi_log_like, y = y)

