#*******************************************************************************
# BASELINE MODEL
#*******************************************************************************

#1. SIMULATE FROM BASE MODEL
#' Simulate an epidemic outbreak
#'
#' Returns daily infection counts from a simulated epidemic outbreak
#'
#' @param num_days Number of days of the epidemic
#' @param shape_gamma Shape of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param scale_gamma Scale of the gamma distribution representing the time-varying infectivity curve of each infected individual
#' @param r0  Model parameter \code{"r0"}. The daily infection count from regular spreading is assumed to follow a poisson distribution with rate \code{r0}*\code{\lambda_t} 
#' @return Total infections; Total daily infection counts of length \code{'num_days'}=
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' tot_daily_infections = SIMULATE_BASELINE_EPIDEMIC(num_days = 50, shape_gamma = 6, scale_gamma = 1, r0 = 1.2)S
#'
#' @export
SIMULATE_EPI_BASELINE  <- function(num_days = 50, r0 = 2.0,
                                 shape_gamma = 6, scale_gamma = 1,
                                 epi_data = c(0,0,0), SIM_DATA = TRUE) {
  
  'Baseline simulation model'
  
  #Initialisation 
  tot_daily_infections = vector('numeric', num_days)
  
  if(SIM_DATA){
    tot_daily_infections[1] = 2
  } else {
    tot_daily_infections[1] = epi_data[1] 
  }
 
  
  #Infectiousness Pressure - Sum of all infected individuals infectivety curves 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  #inf_count = 0
  for (t in 2:num_days) {
    
    #Total rate
    tot_rate = r0*sum(tot_daily_infections[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infections & their probablilty of infection along the gamma dist at that point in time
    tot_daily_infections[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  tot_daily_infections
}

#*
#' Log likelihood of the baseline epidemic model
#'
#' Returns the log likelihood of the baseline epidemic model for given \code{"epidemic_data"} and reproduction number \code{"r0"}
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param r0 Reproduction number r0
#' @return Log likelihood of the baseline model
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#' r0 = 1.2
#' log_likelihood_baseline = LOG_LIKE_BASELINE(epidemic_data, r0)
#'
#' @export
LOG_LIKE_BASELINE <- function(epidemic_data, r0){
  
  #Params
  num_days = length(epidemic_data)
  
  #Infectivity of each individual  - shape and scale of gamma distribution representing the infectivity of each individual 
  shape_gamma = 6; scale_gamma = 1
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  log_likelihood = 0; total_rate = 0
  
  for (t in 2:num_days) {
    
    lambdaX = r0*sum(epidemic_data[1:t-1]*rev(prob_infect[1:t-1]))
    #log_likelihood1 = log_likelihood + epidemic_data[t]*log(lambdaX) - lambdaX - lfactorial(epidemic_data[t]) #Need to include normalizing constant 
    log_likelihood = log_likelihood + dpois(epidemic_data[t], lambda = lambdaX, log = TRUE)

  }
  
  log_likelihood
}

#******************
# LOG LIKE POINT WISE
LOG_LIKE_BASELINE_POINTWISE <- function(epidemic_data, r0){
  
  #Params
  num_days = length(epidemic_data)
  log_like_vec = length(epidemic_data-1)
  
  #Infectivity of each individual  - shape and scale of gamma distribution representing the infectivity of each individual 
  shape_gamma = 6; scale_gamma = 1
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  total_rate = 0
  
  for (t in 2:num_days) {
    
    lambdaX = r0*sum(epidemic_data[1:t-1]*rev(prob_infect[1:t-1]))
    #log_likelihood1 = log_likelihood + epidemic_data[t]*log(lambdaX) - lambdaX - lfactorial(epidemic_data[t]) #Need to include normalizing constant 
    log_like_vec[t-1] = dpois(epidemic_data[t], lambda = lambdaX, log = TRUE)
    
  }
  
  return(log_like_vec)
}

#****************
#* SET PRIOR SSE
#* **************
SET_BASELINE_PRIOR <- function(r0, r0_dash, PRIORS_USED){
  
  #PRIOR
  list_priors = GET_LIST_PRIORS_BASELINE() 
  
  if (PRIORS_USED$BASELINE$r0$EXP) {
    
    prior = dexp(r0_dash, rate = list_priors$exp[1], log = TRUE) -
      dexp(r0, rate = list_priors$exp[1], log = TRUE) 
    
  } else if (PRIORS_USED$BASELINE$r0$GAMMA){
    #print('Gamma prior used r0 sse')
    
    prior = dgamma(r0_dash, shape = list_priors$gamma[1], scale = list_priors$gamma[2], log = TRUE) -
      dgamma(r0, shape = list_priors$gamma[1], scale = list_priors$gamma[2], log = TRUE) 
  }
  
  else if (PRIORS_USED$BASELINE$r0$UNIF) {
    
    prior = dunif(r0_dash, min = list_priors$unif[1], max = list_priors$unif[2], log = TRUE) - 
      dunif(r0, min = list_priors$unif[1], max = list_priors$unif[2], log = TRUE)

  }
  
  return(prior)
}

#************************************************************************
# BASELINE MCMC
#************************************************************************
#'
#' MCMC adaptive algorithm for a baseline epidemic model
#'
#' MCMC algorithm with adaptation for obtaining samples of the reproduction number r0 of an epidemic model
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of MCMC specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the MCMC sampler.
#'   \item \code{"r0_start"} - Model parameter starting points; where the MCMC algorithm begins sampling the reproduction number r0 from.
#'   \item \code{"r0_prior_exp"} - Prior on r0. Default prior is an exponential distribution with rate 1.
#'   \item \code{"target_accept_rate"} - target acceptance rate. Used in the adaptive algorithm
#'   \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept.
#' }
#' @param FLAGS_LIST A list of Boolean variables for switching on/off certain functionality
#' \itemize{
#'   \item \code{"ADAPTIVE"} - Adaptive MCMC algorithm implemented if TRUE
#'   \item \code{"PRIOR"}  - Apply prior distributions to model parameters
#'   \item \code{"THIN"}  - Return a thinned MCMC sample if TRUE, reduced by a factor of \code{"thinning_factor"}
#' }
#' @return MCMC samples, acceptance rates etc.
#' \itemize{
#'   \item \code{"r0_vec"} - A vector containing MCMC samples of r0, the repoduction number of the epidemic.
#'   \item \code{"log_like_vec"}  - A vector containing the log likelihood at each of the MCMC iterations, of size of \code{"mcmc_vec_size"}.
#'   \item \code{"sigma_vec"} - A vector containing sigma adapted at each iteration of the MCMC algorithm, of size of \code{"mcmc_vec_size"} if \code{"FLAGS_LIST$ADAPTIVE"} is TRUE. Otherwise sigma_vec is a constant value.
#'   \item \code{"accept_rate"}  - The MCMC acceptance rate (percentage) for r0. 
#'   }
#'@export
#'
#'@author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' mcmc_inputs = list(n_mcmc = 500000, r0_start = 1.1, r0_prior_exp = c(1, 0), target_accept_rate = 0.4, thinning_factor = 10)
#'
#' #START MCMC
#' mcmc_baseline_output = SSE_MCMC_ADAPTIVE(epidemic_data, mcmc_inputs)

GET_DIC <- function(log_lik_matrix) {
  # log_lik_matrix is a matrix of log-likelihood values
  # Rows correspond to MCMC samples, columns correspond to data points
  
  mean_log_lik <- mean(rowSums(log_lik_matrix))
  
  deviance_mean <- -2 * mean_log_lik
  
  deviance_mcmc <- -2 * rowSums(log_lik_matrix)
  
  p_dic <- 0.5 * var(deviance_mcmc)
  
  dic <- deviance_mean + p_dic
  
  #return(list(DIC = dic, deviance_mean = deviance_mean, p_dic = p_dic))
  
  return(dic)
}


#' @export
MCMC_INFER_BASELINE <- function(epidemic_data, n_mcmc, 
                                PRIORS_USED = GET_PRIORS_USED(), #PRIORS = list(EXP = TRUE, UNIF = FALSE, GAMMA = FALSE),
                                FLAGS_LIST = list(ADAPTIVE = TRUE, 
                                                  THIN = TRUE, BURN_IN = TRUE, COMPUTE_WAIC = TRUE),
                                mcmc_inputs = list(target_accept_rate = 0.4, thinning_factor = 10, #r0_start = 1.0, #1.2 changed 24/12/23
                                                      burn_in_pc = 0.2), B = 500) {
  
  'Returns MCMC samples of the reproduction number of the data
  and acceptance rates'
  
  'Prior
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a'
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC INITIAL POINTS
  r0_start = GET_R0_INITIAL_MCMC(epidemic_data)

  #THINNING FACTOR
  print(paste0('num mcmc iters = ', n_mcmc))
  i_thin = 1
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #BURN-IN
  if(FLAGS_LIST$BURN_IN){
    burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; print(paste0('N burn-in = ', burn_in_start))
    mcmc_vec_size =  mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size; ; print(paste0('Post burn-in mcmc vec size  = ', mcmc_vec_size))
  }
  
  #MCMC VECTORS - INITIALISE
  r0_vec <- vector('numeric', mcmc_vec_size); 
  log_like_vec <- vector('numeric', mcmc_vec_size)
  loglike_pointwise_matrix = matrix(nrow = mcmc_vec_size, ncol = length(epidemic_data)-1)
  r0_vec[1] <- r0_start
  r0 = r0_vec[1]; 
  
  log_like_vec[1] <- LOG_LIKE_BASELINE(epidemic_data, r0_vec[1])
  log_likelihood = log_like_vec[1]
  count_accept = 0
  
  #SIGMA - INITIALISE
  sigma_i =  0.25 #*r0_start
  
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma_vec <- vector('numeric', mcmc_vec_size); sigma_vec[1] =  sigma_i;
    #Other adaptive parameters
    delta = 1/(mcmc_inputs$target_accept_rate*(1-mcmc_inputs$target_accept_rate))
  } else {
    sigma_vec = sigma_i
  }
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    #****************************************************** s
    #r0 proposal
    r0_dash <- r0 + rnorm(1, sd = sigma_i)
    
    if(r0_dash < 0){
      r0_dash = abs(r0_dash)
    }
    
    loglike_new = LOG_LIKE_BASELINE(epidemic_data, r0_dash)
    log_accept_ratio = loglike_new - log_likelihood 
    
    #PRIORS
    log_accept_ratio = log_accept_ratio + SET_BASELINE_PRIOR(r0, r0_dash, PRIORS_USED)
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      r0 <- r0_dash
      count_accept = count_accept + 1
      log_likelihood = loglike_new
    }
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio)) #Acceptance PROB = MIN(1, EXP(ACCPET_PROB))
      #sigma_i =  sigma_i*exp(delta/max(1,1+i-B)*(accept_prob - mcmc_inputs$target_accept_rate))
      sigma_i =  sigma_i*exp(delta/(1+i)*(accept_prob - mcmc_inputs$target_accept_rate))
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {

      r0_vec[i_thin] <- r0
      log_like_vec[i_thin] <- log_likelihood
      
      if (FLAGS_LIST$ADAPTIVE){
        sigma_vec[i_thin] = sigma_i
      }
      
      if (FLAGS_LIST$COMPUTE_WAIC){
        ll_vec = LOG_LIKE_BASELINE_POINTWISE(epidemic_data, r0)
        loglike_pointwise_matrix[i_thin,] = ll_vec
      }
      
      i_thin = i_thin + 1
    }
  }
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  print(paste0('accept_rate = ', accept_rate))
  
  #COMPUTE WAIC, DIC
  waic_result = WAIC(loglike_pointwise_matrix)$WAIC
  dic_result = GET_DIC(loglike_pointwise_matrix)
  
  #Return a, acceptance rate
  return(list(r0_vec = r0_vec, log_like_vec = log_like_vec, sigma_vec = sigma_vec,
              waic_result = waic_result, dic_result = dic_result,
              accept_rate = accept_rate, r0_start = r0_start))
}
