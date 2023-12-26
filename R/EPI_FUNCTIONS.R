#EPIDEMIC FUNCTIONS

#R0 - INITIAL CONDITIONS
GET_R0_INITIAL_MCMC <-function(epi_data){
  
  'INITIALSE R0 DEPENDING ON DATA'
  
  if(sum(epi_data)<=100){
    r0_start = runif(1, 0.9, 1.2)
  } else if( sum(epi_data)> 100 && sum(epi_data)<= 1000){
    r0_start = runif(1, 1, 2)
  }  else if( sum(epi_data)>= 1000 && sum(epi_data)<= 10000){
    r0_start = runif(1, 1, 3)
  } else if (sum(epi_data) > 10000){
    r0_start = runif(1, 2, 4)
  }
  
  return(r0_start)
}

#LAMBDA FUNCTION -SUM OF INFECTIVITY
#' @export
get_lambda <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  #Parameters
  num_days = length(epidemic_data)
  lambda_vec = vector("numeric", num_days)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Lambda -> days of the infection
  for (t in 1:num_days){
    print(paste0('t = ',t))
    lambda_vec[t] = sum(epidemic_data[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    print(lambda_vec[t])
  }
  
  return(lambda_vec)
}

get_infectious_curve <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  '#Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between '
  #Parameters
  num_days = length(epidemic_data)
  infect_vec = vector("numeric", num_days)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  #print(paste0('prob_infect = ', prob_infect))
  #Lambda -> days of the infection
  for (t in 1:num_days){
    print(paste0('t = ',t))
    print(rev(prob_infect[1:(t-1)]))
    infect_vec[t] = rev(prob_infect[1:(t-1)])
    print(infect_vec[t])
  }
  
  return(infect_vec)
}
