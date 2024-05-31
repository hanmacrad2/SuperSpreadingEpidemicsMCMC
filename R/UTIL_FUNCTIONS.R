#********************
#UTIL FUNCTIONS

#TIME FUNCTIONS
#' @export
GET_CURRENT_TIME_STAMP <- function(){
  
  current_time <- Sys.time()
  time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  
  return(time_string)
}
#' @export
GET_FOLDER_TIME_STAMP <- function(folder_type){
  
  current_time <- Sys.time()
  time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  CURRENT_FOLDER <- paste0(folder_type, '_', time_string, '/')
  print(CURRENT_FOLDER)
  
  return(CURRENT_FOLDER)
}

GET_FOLDER_TIME_STAMPV0 <- function(folder_type, array_index){
  
  current_time <- Sys.time()
  time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  CURRENT_FOLDER <- paste0(folder_type, '_', array_index, '_', time_string, '/')
  print(CURRENT_FOLDER)
  
  return(CURRENT_FOLDER)
}

#' @export #Time elapsed
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
#' @export
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
# LOG_SUM_EXP <- function(vectorX){
#   
#   #REMOVE NA VALUES
#   vectorX = na.omit(vectorX)
#   
#   max_val = max(vectorX)
#   
#   out = max_val + log(sum(exp(vectorX - max_val)))
#   
#   return(out)
# }

#Example
# ll = c(-1000, -1000, -1000)
# lse = LOG_SUM_EXP(ll)
# ll_norm = exp(ll - lse) 
# ll_norm

#************
#* PLOTTING FUNCTIONS
#' @export
PLOT_CUM_MEAN <- function(mcmc_chain, titleX = 'MCMC chain', ylabX = 'MCMC chain') {
  
  #MEAN
  #titleX = paste0(titleX, ' mean'); ylabX = paste0(ylabX, ' mean'); 
  mcmc_mean = cumsum(mcmc_chain)/seq_along(mcmc_chain)
  plot(seq_along(mcmc_mean), mcmc_mean,
       xlab = 'Time', ylab = ylabX,
       main = titleX, #bquote(bold(R[0] ~ "MCMC mean, Start:" ~ .(mcmc_specs$mod_start_points$m1))),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}
#CUM MEAN MCMC
PLOT_CUM_MEAN_MCMC <- function(mcmc_chain, titleX = 'MCMC chain', ylabX = 'MCMC chain',
                               ylims = c(0, 5)) {
  
  #MEAN
  #titleX = paste0(titleX, ' mean'); ylabX = paste0(ylabX, ' mean'); 
  mcmc_mean = cumsum(mcmc_chain)/seq_along(mcmc_chain)
  plot(seq_along(mcmc_mean), mcmc_mean,
       xlab = 'Time', ylab = ylabX,
       ylim = ylims,
       main = titleX, #bquote(bold(R[0] ~ "MCMC mean, Start:" ~ .(mcmc_specs$mod_start_points$m1))),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}

#EPIDEMIC FUNCTIONS
#' @export
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

#GET INFECTIOUS GAMMA CURVE
GET_INFECT_GAMMA_CURVE <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  #Parameters
  num_days = length(epidemic_data)
  infect_curve_ga = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  return(infect_curve_ga)
}

#infect_curve_ga = get_infect_ga_curve(epi_data)

#GET INFECTIVITY
#' @export
get_infectious_curve <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  '#Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between '
  #Parameters
  num_days = length(epidemic_data)
  infect_vec = vector("numeric", num_days)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Lambda -> days of the infection
  for (t in 1:num_days){
    infect_vec[t] = rev(prob_infect[1:(t-1)])
  }
  
  return(infect_vec)
}

#MCMC FUNCTIONS

#**********
#CREDIBLE INTERVALS

#MCMC COLUMNS
#' @export
get_lower_ci <- function(mcmc_col, level = 0.95){
  
  #print(CredibleInterval(mcmc_col, level = 0.95)[[2]])
  mcmc_col[is.na(mcmc_col)] <- 0
  lower_interval = CredibleInterval(mcmc_col, level = level)[[2]]
  
  return(lower_interval)
}

#' @export
get_upper_ci <- function(mcmc_col, level = 0.95){
  
  mcmc_col[is.na(mcmc_col)] <- 0
  #print(CredibleInterval(mcmc_col, level = 0.95)[[3]])
  upper_interval = CredibleInterval(mcmc_col, level = level)[[3]]
  
  return(upper_interval)
}

# 'GET MEAN'
get_mean <- function(data_col, rounding_number = 2){

  mean_col = round(mean(data_col), rounding_number)
  
  print(paste0('Mean: ', mean_col))
  
  return(mean_col)
}

# 'GET MEAN AND CIS'
get_mean_cis <- function(data_col, rounding_number = 2){
  
  lower_ci = round(get_lower_ci(data_col), rounding_number)
  upper_ci = round(get_upper_ci(data_col), rounding_number)
  mean_col = round(mean(data_col), rounding_number)
  
  print(paste0('Mean: ', mean_col, ', CI: [', lower_ci, ', ', upper_ci, ' ]'))
}

get_cis_and_mean <- function(data_vec, rounding_number = 2){
  
  lower_ci = round(get_lower_ci(data_vec), rounding_number)
  upper_ci = round(get_upper_ci(data_vec), rounding_number)
  mean_vec = round(mean(data_vec), rounding_number)
  
  print(paste0('Mean: ', mean_vec, ', CI: [', lower_ci, ', ', upper_ci, ' ]'))
}

#' @export
get_ci_matrix <- function(mcmc_matrix, level = 0.95){
  
  #Storage
  num_cols = dim(mcmc_matrix)[2]
  vec_lower = vector("numeric", length = num_cols)
  vec_upper = vector("numeric", length = num_cols)
  
  for (i in seq(1, num_cols)){
    
    #print(mcmc_matrix[, i])
    
    vec_lower[i] = get_lower_ci(mcmc_matrix[, i], level = level)
    vec_upper[i] = get_upper_ci(mcmc_matrix[, i], level = level)
    
  }
  
  return(list(vec_lower = vec_lower, vec_upper = vec_upper))
}

colMedians <- function(matrix_values) {
  apply(matrix_values, 2, median)
}

#TIME 
#CURRENT DATE TIME FOLDER
GET_FOLDER_TIME_STAMP <- function(folder_type){
  
  current_time <- Sys.time()
  time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  CURRENT_FOLDER <- paste0(folder_type, time_string)
  print(CURRENT_FOLDER)
  
  return(CURRENT_FOLDER)
}


# current_time <- Sys.time()
# 
# # Format the current time into a string
# time_string <- format(current_time, "%Y-%m-%d_%H-%M-%S")
# 
# # Create the folder name using the formatted time string
# CURRENT_FOLDER <- paste0('data_base_', time_string)
# 
# # Print the folder name
# print(CURRENT_FOLDER)

#****************
# FILES
#****************

#SORTED FILES
#numeric_part <- as.numeric(sub("^(\\d+).*", "\\1", basename(filenames)))

# Sort filenames based on numeric part
#sorted_filenames <- filenames[order(numeric_part)]
