#' Simulate an epidemic outbreak
#'
#' Returns daily infection counts from a simulated epidemic outbreak
#'
#' @param num_days Number of days of the epidemic
#' @param shape_gamma Shape of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param scale_gamma Scale of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param R0  Model parameter \code{"R0"}. The daily infection count from regular spreading is assumed to follow a poisson distribution with rate \code{R0}*\code{\lambda_t} 
#' @return Total infections; Total daily infection counts of length \code{'num_days'}
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' tot_daily_infections = SIMULATE_BASELINE_EPIDEMIC(num_days = 50, shape_gamma = 6, scale_gamma = 1, R0 = 1.2)S
#'
SIMULATE_BASELINE_EPIDEMIC = function(R0, num_days = 50, shape_gamma = 6, scale_gamma = 1) {
  
  'Baseline simulation model'
  
  #Initialisation 
  tot_daily_infections = vector('numeric', num_days)
  tot_daily_infections[1] = 2
  
  #Infectiousness Pressure - Sum of all infected individuals infectivety curves 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Total rate
    tot_rate = R0*sum(tot_daily_infections[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infections & their probablilty of infection along the gamma dist at that point in time
    tot_daily_infections[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  tot_daily_infections
}

#' Simulate an epidemic with super-spreading events \code{'SSE'}
#'
#' Returns daily infection counts from a simulated epidemic outbreak
#'
#' @param num_days Number of days of the epidemic
#' @param shape_gamma Shape of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param scale_gamma Scale of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param alpha_X  Model parameter \code{"alpha"}. The daily infection count from regular spreading is assumed to follow a poisson distribution with rate \code{alpha}*\code{\lambda_t} 
#' @param beta_X SSE model parameter \code{"beta"}. 
#' @param gamma_X SSE model parameter \code{"gamma"}
#' @return Total infections; Total daily infection counts of length \code{'num_days'}
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' tot_daily_sse_infections = SIMULATE_SS_EVENTS(num_days = 50, shape_gamma = 6, scale_gamma = 1,
#'  alphaX = 0.8, betaX = 0.02, gammaX = 20)
#'
SIMULATE_SS_EVENTS_NB = function(alphaX, betaX, gammaX, num_days = 50, shape_gamma = 6, scale_gamma = 1) {
  
  #Initialise parameters 
  tot_daily_sse_infections = vector('numeric', num_days)
  nsse_infections = vector('numeric', num_days); sse_infections = vector('numeric', num_days)
  tot_daily_sse_infections[1] = 2; nsse_infections[1] = 2; sse_infections[1] = 0
  
  #Infectiousness Pressure - Sum of all infected individuals infectivety curves 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infections
    lambda_t = sum(tot_daily_sse_infections[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Eg. Today is day 10. Infected on day 1. Current Infectiousness is t - day_infected 
    tot_rate = alphaX*lambda_t #Product of infections & their probablilty of infection along the gamma dist at that point in time
    nsse_infections[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data. *Rate of poisson dist is sum of alpha*lambta_t for each individual
    
    #Super-spreading events
    sse_infections[t] = rnbinom(1, betaX*lambda_t, 1/(1 + gammaX)) #z_t: Total infections due to super-spreading event - num of events x num individuals
    
    #Total infections
    tot_daily_sse_infections[t] = nsse_infections[t] + sse_infections[t]
  }
  
  tot_daily_sse_infections
}

#' Simulate an epidemic with super-spreading individuals \code{'SSI'}
#'
#' Returns daily infection counts from a simulated epidemic outbreak with super-spreading individuals
#'
#' @param num_days Number of days of the epidemic
#' @param shape_gamma Shape of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param scale_gamma Scale of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param alpha_X  Model parameter \code{"alpha"}. The daily infection count from regular spreading is assumed to follow a poisson distribution with rate \code{alpha}*\code{\lambda_t} 
#' @param beta_X SSE model parameter \code{"beta"}. 
#' @param gamma_X SSE model parameter \code{"gamma"}
#' @return Total infections; Total daily infection counts of length \code{'num_days'}
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' tot_daily_ssi_infections = SIMULATE_SS_INDIVIDUALS(num_days = 50, shape_gamma = 6, scale_gamma = 1,
#'  alphaX = 0.8, betaX = 0.02, gammaX = 20)
#'
SIMULATE_SS_INDIVIDUALS = function(num_days, shape_gamma, scale_gamma, aX, bX, cX) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Initialise infections
  tot_daily_ssi_infections = vector('numeric', num_days)
  nssi_infections = vector('numeric', num_days); ssi_infections = vector('numeric', num_days)
  tot_daily_ssi_infections[1] = 3; nssi_infections[1] = 2; ssi_infections[1] = 1 
  
  #Infectiousness Pressure - Sum of all infected individuals infectivety curves 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infections
    lambda_t = sum((tot_daily_ssi_infections[1:(t-1)] + cX*ssi_infections[1:(t-1)])*rev(prob_infect[1:(t-1)])) #Product of infections & their probablilty of infection along the gamma dist at that point in time
    nssi_infections[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ssi_infections[t] = rpois(1, bX*lambda_t)
    tot_daily_ssi_infections[t] = nssi_infections[t] + ssi_infections[t]
  }
  
  return(list(nssi_infections, ssi_infections))
}

#************
SIMULATE_SSE_BRANCHING = function(shape_gamma, scale_gamma, alphaX, betaX, gammaX) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nsse_infecteds = vector('numeric', num_days)
  sse_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 2
  nsse_infecteds[1] = 2
  sse_infecteds[1] = 0
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(total_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, betaX*lambda_t) #Number of super-spreading events (beta)
    sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}