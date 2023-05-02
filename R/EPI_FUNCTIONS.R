#EPIDEMIC FUNCTIONS

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
    lambda_vec[t] = sum(epidemic_data[1:(t-1)]*rev(prob_infect[1:(t-1)]))
  }
  
  return(lambda_vec)
}