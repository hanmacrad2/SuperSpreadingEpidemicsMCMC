#*******************************************************************************
# BASELINE MODEL
#*******************************************************************************
#*
#' Log likelihood of the baseline epidemic model
#'
#' Returns the log likelihood of the baseline epidemic model for given \code{"epidemic_data"} and reproduction number \code{"R0"}
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param R0 Reproduction number R0
#' @return Log likelihood of the baseline model
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#' R0 = 1.2
#' log_likelihood_baseline = LOG_LIKE_BASELINE(epidemic_data, R0)
#'
LOG_LIKE_BASELINE <- function(epidemic_data, R0){
  
  #Params
  num_days = length(epidemic_data)
  
  #Infectivity of each individual  - shape and scale of gamma distribution representing the infectivity of each individual 
  shape_gamma = 6; scale_gamma = 1
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  log_likelihood = 0
  
  for (t in 2:num_days) {
    
    lambda = R0*sum(epidemic_data[1:t-1]*rev(prob_infect[1:t-1]))
    log_likelihood = log_likelihood + epidemic_data[t]*log(lambda) - lambda - lfactorial(epidemic_data[t]) #Need to include normalizing constant 
    
  }
  
  log_likelihood
}

#************************************************************************
# BASELINE MCMC
#************************************************************************
#'
#' MCMC adaptive algorithm for a baseline epidemic model
#'
#' MCMC algorithm with adaptation for obtaining samples of the reproduction number R0 of an epidemic model
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of MCMC specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the MCMC sampler.
#'   \item \code{"r0_start"} - Model parameter starting points; where the MCMC algorithm begins sampling the reproduction number R0 from.
#'   \item \code{"r0_prior_exp"} - Prior on R0. Default prior is an exponential distribution with rate 1.
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
#'   \item \code{"r0_vec"} - A vector containing MCMC samples of R0, the repoduction number of the epidemic.
#'   \item \code{"log_like_vec"}  - A vector containing the log likelihood at each of the MCMC iterations, of size of \code{"mcmc_vec_size"}.
#'   \item \code{"sigma_vec"} - A vector containing sigma adapted at each iteration of the MCMC algorithm, of size of \code{"mcmc_vec_size"} if \code{"FLAGS_LIST$ADAPTIVE"} is TRUE. Otherwise sigma_vec is a constant value.
#'   \item \code{"accept_rate"}  - The MCMC acceptance rate (percentage) for R0. 
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
#'
#'
BASELINE_MCMC_ADAPTIVE <- function(epidemic_data,
                                   mcmc_inputs = list(n_mcmc = 1000, r0_start = 1.1, r0_prior_exp = c(1, 0),
                                                      target_accept_rate = 0.4, thinning_factor = 10), 
                                   FLAGS_LIST = list(ADAPTIVE = TRUE, PRIOR = TRUE, THIN = TRUE)) {
  
  'Returns MCMC samples of the reproduction number of the data
  and acceptance rates'
  
  'Prior
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a'
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  n_mcmc = mcmc_inputs$n_mcmc
  print(paste0('num mcmc iters = ', n_mcmc))
  count_accept = 0
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #MCMC VECTORS - INITIALISE
  r0_vec <- vector('numeric', mcmc_vec_size); log_like_vec <- vector('numeric', mcmc_vec_size)
  r0_vec[1] <- mcmc_inputs$r0_start
  log_like_vec[1] <- LOG_LIKE_BASELINE(epidemic_data, r0_vec[1])
  #Running parameters
  r0 = r0_vec[1]; log_likelihood = log_like_vec[1]
  
  #SIGMA - INITIALISE
  sigmaX =  0.5*mcmc_inputs$r0_start
  
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma_vec <- vector('numeric', mcmc_vec_size); sigma_vec[1] =  sigmaX;
    #Other adaptive parameters
    delta = 1/(mcmc_inputs$target_accept_rate*(1-mcmc_inputs$target_accept_rate))
  } else {
    sigma_vec = sigmaX
  }
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    #****************************************************** s
    #R0 proposal
    r0_dash <- r0 + rnorm(1, sd = sigmaX)
    
    if(r0_dash < 0){
      r0_dash = abs(r0_dash)
    }
    
    loglike_new = LOG_LIKE_BASELINE(epidemic_data, r0_dash)
    log_accept_ratio = loglike_new - log_likelihood 
    
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_ratio = log_accept_ratio - r0_dash + r0 #Acceptance RATIO. Acceptance PROB = MIN(1, EXP(ACCPET_PROB))
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      r0 <- r0_dash
      count_accept = count_accept + 1
      log_likelihood = loglike_new
    }
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio)) #Acceptance PROB = MIN(1, EXP(ACCPET_PROB))
      sigmaX =  sigmaX*exp(delta/(1+i)*(accept_prob - mcmc_inputs$target_accept_rate))
    }
    
    #POPULATE MCMC VECTOR (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      i_thin = i/thinning_factor
      r0_vec[i_thin] <- r0
      log_like_vec[i_thin] <- log_likelihood
      
      if (FLAGS_LIST$ADAPTIVE){
        sigma_vec[i_thin] = sigmaX
      }
    }
  }
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  print(paste0('accept_rate = ', accept_rate))
  
  #Return a, acceptance rate
  return(list(r0_vec = r0_vec, log_like_vec = log_like_vec, sigma_vec = sigma_vec,
              accept_rate = accept_rate))
}
