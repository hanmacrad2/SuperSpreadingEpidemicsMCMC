#' Simulation function for the Super-Spreading Individuals (SSI) epidemic model
#'
#' Returns a time series of data of the super-spreading individuals (SSI)epidemic model for given epidemic \code{"data"} and set of model parameters \code{"a, b"} and \code{"c"} and
#'
#' @export
#' 
SIMULATE_EPI_SSIB = function(num_days = 50, r0 = 2.0, alpha = 0.5, b = 10,
                                  shape_gamma = 6, scale_gamma = 1) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Params
  c = (r0*(1 - alpha))/b #r0 = a_prop*r0 + b*c
  #alpha = alpha_prop*r0 
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 3
  nss_infecteds[1] = 2
  ss_infecteds[1] = 1 
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((nss_infecteds[1:(t-1)] + b*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, alpha*r0*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, c*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  total_infecteds
}
#' Log likelihood of the Super-Spreading Individuals (SSI) epidemic model
#'
#' Returns the Log likelihood of the super-spreading individuals (SSI)epidemic model for given epidemic \code{"data"} and set of model parameters \code{"a, b"} and \code{"c"} and
#'
#' @param data data from the epidemic, namely daily infection counts
#' @param r0 SSI model parameter \code{"r0"}
#' @param alpha SSI model parameter \code{"alpha"}
#' @param c SSI model parameter \code{"c"}
#' @return Log likelihood of the SSI model
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' log_like = LOG_LIKE_SSIB(epidemic_data, 0.8, 0.02, 20)
#'
#' @export

LOG_LIKE_SSIB <- function(epidemic_data, r0, alpha, b, 
                          shape_gamma = 6, scale_gamma = 1, temp = 1.2){
  
  #A = PROPORTION OF r0
  c = (r0*(1 - alpha))/b #r0 = a_prop*r0 + b*c
  a = alpha*r0 
  
  non_ss = epidemic_data[[1]]; ss = epidemic_data[[2]]
  
  #Params
  num_days = length(non_ss)
  loglike = 0

  #INFECTIOUSNESS  - Difference of two GAMMA distributions. Discretized
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) { 
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS
    lambda_t = sum((non_ss[1:(t-1)] + b*ss[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    loglike = loglike + dpois(non_ss[t], a*lambda_t, log = TRUE) +
      dpois(ss[t], c*lambda_t, log = TRUE) 
  }
  
  return(loglike)
}

#PRIOR
SET_SSIB_PRIOR <- function(param, param_dash, PRIORS_USED,
                           alpha_flag = FALSE, r0_flag = FALSE, b_flag = FALSE){
  
  #PRIORS
  list_priors = GET_LIST_PRIORS_SSIB(); 
  
  if(alpha_flag){
    
    #BETA PRIOR ON a
    if (PRIORS_USED$SSIB$alpha$BETA) {
      shape1 = list_priors$alpha[1]
      shape2 = list_priors$alpha[2]
      p = dbeta(param_dash, shape1, shape2, log = TRUE) -
        dbeta(param, shape1, shape2, log = TRUE) 
    }
    
  } else if (r0_flag) {
    
    if (PRIORS_USED$SSIB$r0$EXP) {
      p = dexp(param_dash, log = TRUE) - dexp(param, log = TRUE) 
    }
    
  } else if (b_flag) {
    
    #GAMMA PRIOR ON c
    if (PRIORS_USED$SSIB$b$GAMMA) {
      shape = list_priors$b[1]
      scale = list_priors$b[2]
      p = dgamma(param_dash-1, shape = shape, scale = scale, log = TRUE) -
        dgamma(param-1, shape = shape, scale = scale, log = TRUE) 
    }
  }

  return(p)  
}


#' 
#' @export
#LOGLIKELIHOOD + DATA AUGMENTATION MODEL EVIDENCE

#' MCMC adaptive algorithm for Super-Spreading Individuals epidemic model
#'
#' MCMC algorithm with Adaptation for obtaining samples from the parameters of a
#' super-spreading individuals (SSI) epidemic model. Includes; Adaptation, Data Augmentation, B-C & A-C transform
#'
#' @param data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of mcmc specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
#'   \item \code{"param_starts"} - Model parameter starting points; where the mcmc algorithm begins sampling from
#'   \item \code{"alpha_star"} - target acceptance rate; used in the adaptive algorithm
#'   \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
#' }
#' @param priors_list A list of prior parameters used
#' \itemize{
#'   \item \code{"a_prior_exp"} - rate of exponential prior on a, default rate of 1
#'   \item \code{"b_prior_exp"} - rate of exponential prior on b, default rate of 1
#'   \item \code{"b_prior_ga"} -  shape and scale of gamma distribution prior on b, defaults Ga(10, 0.02)
#'   \item \code{"c_prior_exp"}  - rate of exponential prior on c, default rate of 0.1
#'   \item \code{"c_prior_ga"} -  shape and scale of gamma distribution prior on c, defaults Ga(10, 1)
#' }
#' @param FLAGS_LIST A list of Boolean variables for switching on/off certain functionality
#' \itemize{
#'   \item \code{"ADAPTIVE"} - Adaptive MCMC algorithm implemented if TRUE
#'   \item \code{"DATA_AUG"} - Data Augmentation implemented as part of the SSI model if TRUE
#'   \item \code{"BCA_TRANSFORM"} - Transform update step proposed as part of MCMC algorithm if TRUE **
#'   \item \code{"PRIOR"}  - Apply prior distributions to model parameters
#'   \item \code{"B_PRIOR_GAMMA"}  - A Gamma prior on b if TRUE, otherwise exponential
#'   \item \code{"C_PRIOR_GAMMA"}  - A Gamma prior on c if TRUE, otherwise exponential
#'   \item \code{"THIN"}  - Return a thinned MCMC sample if TRUE, reduced by a factor of \code{"thinning_factor"}
#' }
#' @return MCMC samples, acceptance rates etc.
#' \itemize{
#'   \item \code{"a_vec"} - A vector containing MCMC samples of the SSI model parameter a. Vector size  mcmc_vec_size = \code{"n_mcmc/thinning_factor"}
#'   \item \code{"b_vec"} - A vector containing MCMC samples of the SSI model parameter b. Vector size \code{"mcmc_vec_size"}
#'   \item \code{"c_vec"} - A vector containing MCMC samples of the SSI model parameter c. Vector size \code{"mcmc_vec_size"}
#'   \item \code{"r0_vec"} - A vector containing samples of the SSI model parameter r0 obained via; r0 = a + b*c. Vector size \code{"mcmc_vec_size"}
#'   \item \code{"log_like_vec"}  - A vector containing the log likelihood at each of the mcmc iterations, of size of \code{"mcmc_vec_size"}
#'   \item \code{"sigma"} - A list of vectors containing the adapted sigma for each of the model parameters at each iteration of the MCMC, if \code{"FLAGS_LIST$ADAPTIVE"} is TRUE. Otherwise sigma is a list of constant sigma values for each of the model parameters
#'   \item \code{"list_accept_rates"}  - A list of the MCMC acceptance rates (%) for each of the model parameters
#'   \item \code{"non_ss"}  - The augmented epidemic data corresponding to the non super-spreaders (non-ss) for each iteration of the mcmc
#'   \item \code{"ss"}  - The augmented epidemic data corresponding to the super-spreaders (ss) for each iteration of the mcmc
#'   }
#' @export
#'
#'@author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#'mcmc_inputs = list(n_mcmc = 500000, param_starts = list(m1 = 0.72, m2 = 0.0038, m3 = 22),
#' alpha_star = 0.4, thinning_factor = 10)
#'
#' #START MCMC
#'mcmc_ssi_output = SSI_MCMC_ADAPTIVE(epidemic_data, mcmc_inputs)
#'
#' @export
#************************************************************************
#1. SSIB MCMC                              (W/ DATA AUGMENTATION)
#************************************************************************
MCMC_INFER_SSIB <- function(epidemic_data, n_mcmc = 40000, PRIORS_USED = GET_PRIORS_USED(),
                            param_starts = list(alpha = 0.5, b = 10, r0_start = 2.0),
                            mcmc_inputs = list(alpha_star = 0.3, 
                                               burn_in_pc = 0.2, thinning_factor = 10),
                            FLAGS_LIST = list(THIN = TRUE, ADAPTIVE = TRUE)){
  
  'Returns MCMC samples of SSI-B model parameters (a, b, c, r0 = a + b*c)
  w/ acceptance rates.
  INCLUDES; ADAPTATION, DATA AUGMENTATION, B-C & A-C transform'
  
  'Priors
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  param_starts$r0_start = r0_start
  
  time = length(epidemic_data)
  print(paste0('num mcmc iters = ', n_mcmc))
  
  #DATA: SUPERSPREADING INTIALISATION
  ss = ifelse(epidemic_data > 1, 1, 0) #Initialising ss to be 1 if epi_data > 1. 0 otherwise
  print('ss: '); print(ss)
  non_ss = epidemic_data - ss 
  print('non_ss: '); print(non_ss)
  data = list(non_ss, ss)
  
  #THINNING FACTOR
  i_thin = 1
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
    } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #BURN_IN: Initial samples are not completely valid; the Markov Chain has not stabilized to the stationary distribution. The burn in samples allow you to discard these initial samples that are not yet at the stationary.
  burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; 
  print(paste0('N burn-in = ', burn_in_start))
  mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size
  mcmc_vec_size = mcmc_vec_size - 1 
  print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
  
  #INITIALISE MCMC VECTORS
  alpha_vec <- vector('numeric', mcmc_vec_size); b_vec <- vector('numeric', mcmc_vec_size)
  c_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);
  
  #INITIALISE MCMC[1]
  alpha_vec[1] <- param_starts$alpha; 
  b_vec[1] <- param_starts$b; r0_vec[1] <- param_starts$r0_start
  #b_vec[1] <- param_starts$m2
  log_like_vec[1] <- LOG_LIKE_SSIB(data, r0_vec[1], alpha_vec[1], b_vec[1])
  
  #INITIALISE RUNNING PARAMS
  r0 = r0_vec[1]
  alpha = alpha_vec[1]; b = b_vec[1];  
  log_like = log_like_vec[1]
  c = ((1 - alpha)*r0)/b; c_vec[1] = c
  print(paste0('loglike: ', log_like))
  
  #SIGMA
  sigma1 =0.4 #param_starts$m1
  sigma2 = 0.5 #0.4*r0
  sigma3 = 0.5*param_starts$b; 
  #sigma4 = 0.85*param_starts$m3; sigma5 = 0.85*param_starts$m3
  
  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){
    
    sigma1_vec <- vector('numeric', mcmc_vec_size); sigma2_vec <- vector('numeric', mcmc_vec_size)
    sigma3_vec <- vector('numeric', mcmc_vec_size)
    sigma1_vec[1] =  sigma1; sigma2_vec[1] =  sigma2; 
    sigma3_vec[1] =  sigma3
    sigma = list()
    #Other adaptive parameters
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
    
  } else {
    sigma = list(sigma1 <- sigma1, sigma2 <- sigma2,
                 sigma3 <- sigma3) #, sigma4 <- sigma4, sigma5 <- sigma5)
  }
  
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0)
                            
  #DATA AUG OUTPUT
  #mat_count_da = matrix(0, mcmc_vec_size, time) #i x t
  non_ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    if (i%%5000 == 0) {
      print(paste0('i = ', i))
    }
    
    #******************************************************
    #alpha
    alpha_dash <- alpha + rnorm(1, sd = sigma1)
    
    #Constraint: 0 < a < 1 (Proporition of non ss)
    while(alpha_dash < 0 || alpha_dash > 1){
      
      if (alpha_dash > 1){  
        alpha_dash = 2 - alpha_dash
      }
      alpha_dash = abs(alpha_dash) 
    }
    
    logl_new = LOG_LIKE_SSIB(data, r0, alpha_dash, b) #RO 
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    
    #PRIOR
    log_accept_ratio = log_accept_ratio + SET_SSIB_PRIOR(alpha, alpha_dash, PRIORS_USED, alpha_flag = TRUE)
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_ratio) {
      alpha <- alpha_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    } 
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma1 =  sigma1*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************************************************ Only if (b > 0) ?
    #r0
    r0_dash <- abs(r0 + rnorm(1, sd = sigma2))

    #loglikelihood
    logl_new = LOG_LIKE_SSIB(data, r0_dash, alpha, b)
    log_accept_ratio = logl_new - log_like
    
    #PRIOR
    log_accept_ratio = log_accept_ratio + SET_SSIB_PRIOR(r0, r0_dash, PRIORS_USED, r0_flag = TRUE)
    #Metropolis Acceptance Step 
    if(log(runif(1)) < log_accept_ratio) {
      r0 <- r0_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    } 
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma2 =  sigma2*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************************************************
    #b
    b_dash <- b + rnorm(1, sd = sigma3)
    
    if(b_dash < 1){
      b_dash = 2 - b_dash #Prior on b: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSIB(data, r0, alpha, b_dash)
    log_accept_ratio = logl_new - log_like
    
    #PRIOR
    log_accept_ratio = log_accept_ratio + SET_SSIB_PRIOR(b, b_dash, PRIORS_USED, b_flag = TRUE)
    
    #Metropolis Acceptance Step
    #if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
    if(log(runif(1)) < log_accept_ratio) {
      b <- b_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma3 =  sigma3*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************
    #DATA AUGMENTATION
    #***********************************
    #FOR EACH S_T
    for(t in 1:time){
      
      data_dash = data
      
      #Propose
      x = data[[1]][t] + data[[2]][t]
      d = round(runif(1, min = 0, max = x/2))
      
      #STOCHASTIC PROPOSAL for s
      if (runif(1) < 0.5) {
        st_dash = data[[2]][t] + d #1 
      } else {
        st_dash = data[[2]][t] - d #1 
      }
      
      #ACCEPTANCE PROBABILITY
      data_dash[[2]][t] = st_dash #s_t = st_dash
      data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t
      
      #CRITERIA FOR S_T & N_T
      if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
        log_accept_ratio = -Inf 
      } else {
        logl_new = LOG_LIKE_SSIB(data_dash, r0, alpha, b) #+include alpha prior 
        log_accept_ratio = logl_new - log_like
      }

      #METROPOLIS ACCEPTANCE STEP
      if(log(runif(1)) < log_accept_ratio) {
        data <- data_dash
        log_like <- logl_new
      }
      
      #Store
      if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
        non_ss[i_thin, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i_thin, t] = data[[2]][t]
      }
    }
    
    #Loglikelihood Check (Passing - no error)
    #if (!(log_like && LOG_LIKE_SSIB(data, a, b, c))) print('loglike doesnt exist')
    #else (log_like!=LOG_LIKE_SSIB(data, a, b, c)) print(paste0('ERROR! logl diff = ', log_like - LOG_LIKE_SSIB(data, a, b, c)))
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin < mcmc_vec_size) {
      alpha_vec[i_thin] <- alpha; r0_vec[i_thin] <- r0
      b_vec[i_thin] <- b; c_vec[i_thin] <- (r0*(1-alpha))/b     #(r0-a)/c #a + b*c
      log_like_vec[i_thin] <- log_like 
      sigma$sigma1[i_thin] = sigma1; sigma$sigma2[i_thin] = sigma2; sigma$sigma3[i_thin] = sigma3
      i_thin = i_thin + 1
    }
  }
  #FINAL
  alpha_vec = alpha_vec[1:length(alpha_vec)-1]; b_vec = b_vec[1:length(b_vec)-1]
  c_vec = c_vec[1:length(c_vec)-1]; r0_vec = r0_vec[1:length(r0_vec)-1]
  log_like_vec = log_like_vec[1:length(log_like_vec)-1]
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1)
  #accept_rate6 = 100*list_accept_counts$count_accept6/((n_mcmc-1)*time) #i x t
  
  #Acceptance rates
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3)
  print(list_accept_rates)
  
  #Return a, acceptance rate
  return(list(alpha_vec = alpha_vec, b_vec = b_vec, c_vec = c_vec, r0_vec = r0_vec,
              log_like_vec = log_like_vec, sigma = sigma,
              list_accept_rates = list_accept_rates,
              data = data, non_ss = non_ss, ss = ss, r0_start = r0_start))
}

