#SIMULATION FUNCTIONS

#***************************************
#BASELINE MODEL
SIMULATE_EPI_BASELINE_REAL  <- function(real_epi_data, num_days = 50, r0 = 2.0,
                                   shape_gamma = 6, scale_gamma = 1,
                                   epi_data = c(0,0,0), REAL_DATA = TRUE) {
  
  'Baseline simulation model'
  
  #Initialisation 
  tot_daily_infections = vector('numeric', num_days)
  
  if(REAL_DATA){
    tot_daily_infections[1] = real_epi_data[1]
  } else {
    tot_daily_infections[1] = 2
  }
  
  #Infectiousness Pressure - Sum of all infected individuals infectivety curves 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Total rate
    tot_rate = r0*sum(tot_daily_infections[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infections & their probablilty of infection along the gamma dist at that point in time
    tot_daily_infections[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  tot_daily_infections
}


#***************************************
#SSE MODEL
SIMULATE_EPI_SSE_REAL <- function(real_epi_data, num_days = 50, r0 = 2.0, k = 0.2, #k = 0.16, 0.8
                             shape_gamma = 6, scale_gamma = 1,
                             epi_data = c(0,0,0), REAL_DATA = TRUE,
                             FLAG_NEGBIN_PARAMATERISATION = list(param_mu = TRUE, param_prob = FALSE)) {
  
  'Simulate an epidemic with Superspreading events
  alpha - r0'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); 
  
  if(REAL_DATA){
    x[1] = real_epi_data[1]
  } else {
    x[1] = 2
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infections (tot_rate = lambda) fix notation
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    
    #NEGATIVE BINOMIAL PARAMETERISATION
    if (FLAG_NEGBIN_PARAMATERISATION$param_mu){
      
      x[t] = rnbinom(1, size = k*lambda_t, mu =  r0*lambda_t) #Neg Bin parameterisation #1
      
    } else if (FLAG_NEGBIN_PARAMATERISATION$param_prob) {
      
      x[t] = rnbinom(1, size = k*lambda_t, prob =  k/(r0 + k)) #Neg Bin parameterisation #2
      
    }
  }
  
  return(x)
}


#*******************************************
#SSI MODEL
SIMULATE_EPI_SSI_REAL <- function(real_epi_data, num_days = 50, r0 = 2.0, k = 0.8,
                             shape_gamma = 6, scale_gamma = 1,
                             epi_data = c(0,0,0), REAL_DATA = TRUE) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); 
  eta_vec = vector('numeric', num_days); 
  
  if(REAL_DATA){
    x[1] = real_epi_data[1]
  } else {
    x[1] = 2 
  }
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Infectiousness Pressure
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #ETA (t-1)
    eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = r0/k) #Draw eta from previous time step
    
    #INFECTIVITY
    infectivity = rev(prob_infect[1:t-1]) 
    
    #POISSON; OFFSPRING DISTRIBUTION
    total_rate = sum(eta_vec[1:t-1]*infectivity) 
    
    x[t] = rpois(1, total_rate)
    
  }
  #result = (list(epidemic_data = x, eta_vec = eta_vec))
  
  return(x) #result
}


#***************************************
#SSEB MODEL
SIMULATE_EPI_SSEB_REAL <- function(real_epi_data, num_days = 50, r0 = 2.0, alpha = 0.5, beta = 10,
                              shape_gamma = 6, scale_gamma = 1,
                              epi_data = c(0,0,0), REAL_DATA = TRUE) {
  
  'Simulate an epidemic with Superspreading events
  alpha - rate of non super-spreading events/days
  Beta = Proportion of superspreading events/days
  beta = increased rate of superspreading event'
  
  #MODEL PARAMS
  gamma = r0*(1 - alpha)/beta #rate of infections = r0_sse/num infections
  #alpha = alpha*r0 #alpha is now the rate
  
  print(paste0('alpha = ', alpha)); print(paste0('beta = ', beta))
  print(paste0('beta = ', beta))
  #Set up
  total_infections = vector('numeric', num_days)
  nsse_infections = vector('numeric', num_days)
  sse_infections = vector('numeric', num_days)
  
  if (REAL_DATA){
    total_infections[1] = real_epi_data[1]
    nsse_infections[1] = real_epi_data[1]
    sse_infections[1] = 0 
  } else {
    total_infections[1] = 2
    nsse_infections[1] = 2
    sse_infections[1] = 0 
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infections (tot_rate = lambda) fix notation
    lambda_t = sum(total_infections[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    tot_rate = alpha*r0*lambda_t #Product of infections & their probablilty of infection along the gamma dist at that point in time
    nsse_infections[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, gamma*lambda_t) #Number of super-spreading events (gamma)
    sse_infections[t] = rpois(1, beta*n_t) #z_t: Total infections due to super-spreading event - num of events x Num individuals
    
    total_infections[t] = nsse_infections[t] + sse_infections[t]
  }
  
  total_infections
}

#************************************
#* SSIB MODEL
SIMULATE_EPI_SSIB_REAL = function(real_epi_data, num_days = 50, r0 = 2.0, a = 0.5, b = 10,
                             data_start = c(3,2,1), REAL_DATA = TRUE,
                             shape_gamma = 6, scale_gamma = 1) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Params
  c = (r0*(1 - a))/b #r0 = a_prop*r0 + b*c
  
  #DATA
  total_infections = vector('numeric', num_days)
  non_ss = vector('numeric', num_days)
  ss = vector('numeric', num_days)
  
  if (REAL_DATA){
    data_point1 = real_epi_data[1]
    total_infections[1] = real_epi_data[1]
    ss[1] = ifelse(data_point1 > 1, max(1, round(0.1*data_point1)), 0) #10% of ssib data if > 1, else 1, else 0
    non_ss[1] =   total_infections[1] - ss[1]
    
  } else {
    total_infections[1] = data_start[1] 
    non_ss[1] =  data_start[2] 
    ss[1] =  data_start[3] 
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infections (tot_rate = lambda) fix notation
    lambda_t = sum((non_ss[1:(t-1)] + b*ss[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infections & their probablilty of infection along the gamma dist at that point in time
    non_ss[t] = rpois(1, a*r0*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss[t] = rpois(1, c*lambda_t)
    total_infections[t] = non_ss[t] + ss[t]
  }
  
  total_infections
}
