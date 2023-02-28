#********************
#UTIL FUNCTIONS

#TIME FUNCTIONS
get_time <- function(start_time, end_time, show = TRUE){
  'Print difference between end & start time'
  
  time_elap = round(end_time - start_time, 2)
  if(show){
    print('Time elapsed:') 
    print(time_elap)
  }
  time_elap
}

#FOLDER CREATION
create_folder <- function(folder_name){
  
  ifelse(!dir.exists(file.path(folder_name)),
         dir.create(file.path(folder_name), recursive = TRUE), FALSE)
}

#'Log-sum-exp of a vector
#'
#' Log-sum-exp trick applied to a vector to overcome the issue of under- or overflow.
#'  Often arise when normalizing vectors of log probabilities
#' 
#' @param vectorX A vector of quantities 
#' @return normalized exponentiated vector
#' @export LOG_SUM_EXP
#'
#' @author Hannah Craddock
#'
#' @examples
#'
#' vectorX = LOG_SUM_EXP(vectorX)
#'

#LOG EXP SUM
LOG_SUM_EXP <- function(vectorX){
  
  max_val = max(vectorX)
  
  out = max_val + log(sum(exp(vectorX - max_val)))
  
  return(out)
}

#Example
# ll = c(-1000, -1000, -1000)
# lse = LOG_SUM_EXP(ll)
# ll_norm = exp(ll - lse) 
# ll_norm

#EPIDEMIC FUNCTIONS
get_lambda <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  '#Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between '
  #Parameters
  num_days = length(epidemic_data)
  lambda_vec = vector("numeric", num_days)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Lambda -> days of the infection
  for (t in 1:num_days){
    lambda_vec[t] = sum(epidemic_data[1:(t-1)]*rev(prob_infect[1:(t-1)]))
  }
  
  return(lambda_vec)
}
#GET INFECTIVITY
get_infectivity <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  '#Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between '
  #Parameters
  num_days = length(epidemic_data)
  infect_vec = vector("numeric", num_days)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Lambda -> days of the infection
  for (t in seq(1, num_days)){
    infect_vec[t] = rev(prob_infect[1:(t-1)])
  }
  
  return(infect_vec)
}

#MCMC FUNCTIONS

#**********
#CREDIBLE INTERVALS

#MCMC COLUMNS
get_lower_ci <- function(mcmc_col){
  
  lower_interval = CredibleInterval(mcmc_col, level = 0.95)[[2]]
  
  return(lower_interval)
}

get_upper_ci <- function(mcmc_col){
  
  upper_interval = CredibleInterval(mcmc_col, level = 0.95)[[3]]
  
  return(upper_interval)
}

get_ci_matrix <- function(mcmc_matrix){
  
  #Storage
  num_cols = dim(mcmc_matrix)[2]
  vec_lower = vector("numeric", length = num_cols)
  vec_upper = vector("numeric", length = num_cols)
  
  for (i in seq(1, num_cols)){
    
    vec_lower[i] = get_lower_ci(mcmc_matrix[, i])
    vec_upper[i] = get_upper_ci(mcmc_matrix[, i])
    
  }
  
  return(list(vec_lower = vec_lower, vec_upper = vec_upper))
}