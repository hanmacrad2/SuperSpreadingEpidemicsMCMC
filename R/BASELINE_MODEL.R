#*******************************************************************************
# BASELINE MODEL

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
    log_likelihood = logl + epidemic_data[t]*log(lambda) - lambda - lfactorial([t]) #Need to include normalizing constant 
    
  }
  log_likelihood
}

#' Log likelihood of the baseline epidemic model
#'
#' Returns the log likelihood of the baseline epidemic model for given \code{"epidemic_data"} and reproduction number \code{"R0"}
#'
#' MCMC adaptive algorithm for the baseline epidemic model
#'
#' MCMC algorithm with Adaptation for obtaining samples from the parameters of a
#' super-spreading events (SSE) epidemic model
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of mcmc specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
#'   \item \code{"r0_start"} - Model parameter starting points; where the mcmc algorithm begins sampling the reproduction number R0 from
#'   \item \code{"alpha_star"} - target acceptance rate; used in the adaptive algorithm
#'   \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
#' }
#************************************************************************
#1. SSE MCMC
#************************************************************************
BASELINE_MCMC_ADAPTIVE <- function(epidemic_data,
                                   mcmc_inputs = list(n_mcmc = 500000, r0_start = 1.1,
                                                      thinning_factor = 10, r0_prior_exp = c(1, 0)), 
                                   FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE, PRIOR = TRUE)) {
  
  'Returns MCMC samples of the reproduction number \code{"R0"} of the data
  and acceptance rates'
  
  'Prior
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a'
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  n_mcmc = mcmc_inputs$n_mcmc; print(paste0('num mcmc iters = ', n_mcmc))
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
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
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
      sigmaX =  sigmaX*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #POPULATE MCMC VECTOR (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      i_thin = i/thinning_factor
      r0_vec[i_thin] <- r0
      log_like_vec[i_thin] <- log_like
      
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

#' MCMC adaptive algorithm for Super-Spreading Events epidemic model
#'
#' MCMC algorithm with Adaptation for obtaining samples from the parameters of a
#' super-spreading events (SSE) epidemic model
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of mcmc specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
#'   \item \code{"mod_start_points"} - Model parameter starting points; where the mcmc algorithm begins sampling from
#'   \item \code{"alpha_star"} - target acceptance rate; used in the adaptive algorithm
#'   \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
#' }
mcmc_r0 <- function(data, n, sigma, burn_in, x0 = 1) {
  
  'Returns MCMC samples of R0. 
  Prior on R_0 (or alpha) = exp(1)'
  
  #Initialisation
  r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  
  #MCMC chain
  for(i in 2:n) {
    r0_dash <- r0_vec[i-1] + rnorm(1, sd = sigma) #, mean = 0, sd = sigma_opt)
    if(r0_dash < 0){
      r0_dash = abs(r0_dash)
    }
    
    #Alpha
    log_alpha = log_like(epidemic_data, r0_dash) - log_like(epidemic_data, r0_vec[i-1]) - r0_dash + r0_vec[i-1] #exponential prior
    #log_alpha = log_like(epidemic_data, Y) - log_like(epidemic_data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    #Should include: Likelihood + prior + propogsal density x2 (Previous time step & Current time step)
    
    if (is.na(log_alpha)){
      print('na value')
      sprintf("r0_dash: %i", r0_dash)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) { 
      r0_vec[i] <- r0_dash
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
    }
  }
  #Final stats
  accept_rate = 100*(count_accept/n)
  print(paste0("Acceptance rate = ",accept_rate))
  
  r0_vec = r0_vec[burn_in:n]
  r0_vec
  
  return(list(r0_vec, accept_rate))
}



#' Log likelihood of the Super-Spreading Individuals (SSI) epidemic model
#'
#' Returns the Log likelihood of the super-spreading individuals (SSI)epidemic model for given \code{"epidemic_data"} and set of model parameters \code{"a, b"} and \code{"c"} and
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param alpha_X SSI model parameter \code{"\alpha"}
#' @param beta_X SSI model parameter \code{"\beta"}
#' @param gamma_X SSI model parameter \code{"\gamma"}
#' @return Log likelihood of the SSE model
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' log_like = LOG_LIKE_SSE_LSE(epidemic_data, 0.8, 0.02, 20)
#'
LOG_LIKE_SSE_LSE <- function(x, alphaX, betaX, gammaX){
  
  #Initialise parameters
  num_days = length(x);  logl = 0
  shape_gamma = 6; scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)),  shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    
    if((x[t] == 0) | is.na(x[t])) { #y_t also equal to zero
      
      logl = logl -(alphaX*lambda_t) -
        (betaX*lambda_t*log(gammaX +1))
      
    } else {
      
      #Terms in inner sum
      inner_sum_vec <- vector('numeric', x[t])
      
      for (y_t in 0:x[t]){ #Sum for all values of y_t up to x_t
        
        #Store inner L(x_i) term in vector position
        inner_sum_vec[y_t + 1] = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
                                    lgamma((x[t] - y_t) + (betaX*lambda_t)) - lgamma(betaX*lambda_t) -
                                    lfactorial(x[t] - y_t) - (betaX*lambda_t*log(gammaX +1)) +
                                    (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
        
      }
      
      #Calculate max element in inner vector, for all y_t for alphagiven t, x[t]
      lx_max = max(inner_sum_vec)
      
      #Calculate lse
      lse = lx_max + log(sum(exp(inner_sum_vec - lx_max) ))
      
      #Add to overall log likelihood
      logl = logl + lse
      
    }
  }
  logl
}

#' MCMC adaptive algorithm for Super-Spreading Events epidemic model
#'
#' MCMC algorithm with Adaptation for obtaining samples from the parameters of a
#' super-spreading events (SSE) epidemic model
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of mcmc specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
#'   \item \code{"mod_start_points"} - Model parameter starting points; where the mcmc algorithm begins sampling from
#'   \item \code{"alpha_star"} - target acceptance rate; used in the adaptive algorithm
#'   \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
#' }
#' @param priors_list A list of prior parameters used
#' \itemize{
#'   \item \code{"alpha_prior_exp"} - rate of exponential prior on alpha, default rate of 1
#'   \item \code{"beta_prior_exp"} - rate of exponential prior on beta, default rate of 1
#'   \item \code{"beta_prior_ga"} -  shape and scale of gamma distribution prior on beta, defaults Ga(10, 0.02)
#'   \item \code{"gamma_prior_exp"}  - rate of exponential prior on gamma, default rate of 0.1
#'   \item \code{"gamma_prior_ga"} -  shape and scale of gamma distribution prior on c, defaults Ga(10, 1)
#' }
#' @param FLAGS_LIST A list of Boolean variables for switching on/off certain functionality
#' \itemize{
#'   \item \code{"ADAPTIVE"} - Adaptive MCMC algorithm implemented if TRUE
#'   \item \code{"ABG_TRANSFORM"} - Transform update step proposed as part of MCMC algorithm if TRUE **
#'   \item \code{"PRIOR"}  - Apply prior distributions to model parameters
#'   \item \code{"BETA_PRIOR_GA"}  - A Gamma prior on beta if TRUE, otherwise exponential
#'   \item \code{"GAMMA_PRIOR_GA"}  - A Gamma prior on gamma if TRUE, otherwise exponential
#'   \item \code{"THIN"}  - Return a thinned MCMC sample if TRUE, reduced by a factor of \code{"thinning_factor"}
#' }
#' @return MCMC samples, acceptance rates etc.
#' \itemize{
#'   \item \code{"alpha_vec"} - A vector containing MCMC samples of the SSI model parameter alpha. Vector size \code{"mcmc_vec_size = n_mcmc/thinning_factor"}
#'   \item \code{"beta_vec"} - A vector containing MCMC samples of the SSI model parameter beta. Vector size \code{"mcmc_vec_size"}
#'   \item \code{"gamma_vec"} - A vector containing MCMC samples of the SSI model parameter gamma. Vector size \code{"mcmc_vec_size"}
#'   \item \code{"r0_vec"} - A vector containing samples of the SSI model parameter r0 obtained via; r0 = alpha + beta*gamma. Vector size \code{"mcmc_vec_size"}
#'   \item \code{"log_like_vec"}  - A vector containing the log likelihood at each of the MCMC iterations, of size of \code{"mcmc_vec_size"}
#'   \item \code{"sigma"} - A list of vectors containing the adapted sigma for each of the model parameters at each iteration of the MCMC, if \code{"FLAGS_LIST$ADAPTIVE"} is TRUE. Otherwise sigma is a list of constant sigma values for each of the model parameters
#'   \item \code{"list_accept_rates"}  - A list of the MCMC acceptance rates (percentage) for each of the model parameters
#'   }
#' @export
#'
#'@author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#'mcmc_inputs = list(n_mcmc = 500000, mod_start_points = list(m1 = 0.8, m2 = 0.05, m3 = 10),
#' alpha_star = 0.4, thinning_factor = 10)
#'
#' #START MCMC
#'mcmc_sse_output = SSE_MCMC_ADAPTIVE(epidemic_data, mcmc_inputs)
#'

#************************************************************************
#1. SSE MCMC
#************************************************************************
SSE_MCMC_ADAPTIVE <- function(epidemic_data,
                              mcmc_inputs = list(n_mcmc = 500000,
                                                 mod_start_points = list(m1 = 0.8, m2 = 0.05, m3 = 10), alpha_star = 0.4,
                                                 thinning_factor = 10),
                              priors_list = list(alpha_prior_exp = c(1, 0), beta_prior_ga = c(10, 2/100), beta_prior_exp = c(0.1,0),
                                                 gamma_prior_ga = c(10, 1), gamma_prior_exp = c(0.1,0)),
                              FLAGS_LIST = list(ADAPTIVE = TRUE, ABG_TRANSFORM = TRUE,
                                                PRIOR = TRUE, BETA_PRIOR_GA = FALSE, GAMMA_PRIOR_GA = FALSE,
                                                THIN = TRUE)) {
  
  'Returns MCMC samples of SSI model parameters (alpha, beta, gamma, r0 = alpha + beta*gamma)
  w/ acceptance rates.
  INCLUDES; ADAPTATION, beta-gamma & alpha-gamma transform'