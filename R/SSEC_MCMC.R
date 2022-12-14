#SSEC model

#SIMULATE
SIMULATE_EPI_SSEC <- function(num_days = 50, R0 = 1.2, k = 0.16,
                              shape_gamma = 6, scale_gamma = 1) {
  
  'Simulate an epidemic with Superspreading events
  alpha - R0'

  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    x[t] = rnbinom(1, size = k, mu =  R0*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  x
}

#************************
#* LOG LIKELIHOOD SSEC
#* ***********************
LOG_LIKE_SSEC <- function(x, lambda_vec, R0, k){
  
  #Params
  num_days = length(x); loglike = 0
  
  for (t in 2:num_days) {

    loglike = loglike + dnbinom(1, size = k, mu =  R0*lambda_vec[t], log = TRUE)
    
  }
  
  return(loglike)
}
