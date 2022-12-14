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
    tot_rate = R0*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    x[t] = rnbinom(1, size = k, mu = tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  x
}

#APPLY
ssec_data = SIMULATE_EPI_SSEC()
plot.ts(ssec_data)
