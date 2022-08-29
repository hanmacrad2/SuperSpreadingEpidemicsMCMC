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

  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'


  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  n_mcmc = mcmc_inputs$n_mcmc;
  print(paste0('num mcmc iters = ', n_mcmc))

  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }

  #INITIALISE MCMC VECTORS
  alpha_vec <- vector('numeric', mcmc_vec_size); beta_vec <- vector('numeric', mcmc_vec_size)
  gamma_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);

  #INITIALISE MCMC[1]
  alpha_vec[1] <- mcmc_inputs$mod_start_points$m1; beta_vec[1] <- mcmc_inputs$mod_start_points$m2
  gamma_vec[1] <- mcmc_inputs$mod_start_points$m3; r0_vec[1] <- alpha_vec[1] + beta_vec[1]*gamma_vec[1]
  log_like_vec[1] <- LOG_LIKE_SSE_LSE(epidemic_data, alpha_vec[1], beta_vec[1], gamma_vec[1])

  #INITIALISE RUNNING PARAMS
  alpha = alpha_vec[1]; beta = beta_vec[1]; gamma = gamma_vec[1]; log_like = log_like_vec[1]

  #SIGMA
  sigma1 =  0.4*mcmc_inputs$mod_start_points$m1;  sigma2 = 0.3*mcmc_inputs$mod_start_points$m2
  sigma3 = 0.4*mcmc_inputs$mod_start_points$m3; sigma4 = 0.85*mcmc_inputs$mod_start_points$m3
  sigma5 = 0.85*mcmc_inputs$mod_start_points$m3

  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){

    #SIGMA
    sigma1_vec <- vector('numeric', mcmc_vec_size); sigma2_vec <- vector('numeric', mcmc_vec_size)
    sigma3_vec <- vector('numeric', mcmc_vec_size); sigma4_vec <- vector('numeric', mcmc_vec_size)
    sigma5_vec <- vector('numeric', mcmc_vec_size);

    #SIGMA; INITIALISE FIRST ELEMENT
    sigma1_vec[1] =  sigma1; sigma2_vec[1] =  sigma2; sigma3_vec[1] =  sigma3
    sigma4_vec[1] =  sigma4; sigma5_vec[1] =  sigma5

    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma = list(sigma1_vec = sigma1_vec, sigma2_vec = sigma2_vec, sigma3_vec = sigma3_vec,
                 sigma4_vec = sigma4_vec, sigma5_vec = sigma5_vec)

    #Other adaptive parameters
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))

  } else {

    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma = list(sigma1 = sigma1, sigma2 = sigma2,
                 sigma3 = sigma3, sigma4 = sigma4,
                 sigma5 = sigma5)
  }

  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0)

  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {

    #****************************************************** s
    #alpha
    alpha_dash <- alpha + rnorm(1, sd = sigma1)

    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }

    #log a
    logl_new = LOG_LIKE_SSE_LSE(epidemic_data, alpha_dash, beta, gamma)
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_ratio = log_accept_ratio - alpha_dash + alpha #*Actually this is the Acceptance RATIO. ACCEPTANCE PROB = MIN(1, EXP(ACCPET_PROB))
    }

    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      alpha <- alpha_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    }

    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma1 =  sigma1*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #************************************************************************ Only if (b > 0){ ?
    #beta
    beta_dash <- beta + rnorm(1, sd = sigma2)
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }

    #loglikelihood
    logl_new = LOG_LIKE_SSE_LSE(epidemic_data, alpha, beta_dash, gamma)
    log_accept_ratio = logl_new - log_like

    #Priors
    if (FLAGS_LIST$BETA_PRIOR_GA){
      log_accept_ratio = log_accept_ratio +
        dgamma(beta_dash, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE) -
        dgamma(beta, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - beta_dash + beta
    }

    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      beta <- beta_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    }

    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma2 =  sigma2*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #************************************************************************
    #gamma
    gamma_dash <- gamma + rnorm(1, sd = sigma3)
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSE_LSE(epidemic_data, alpha, beta, gamma_dash)
    log_accept_ratio = logl_new - log_like

    #Priors
    if(FLAGS_LIST$GAMMA_PRIOR_GA){
      log_accept_ratio = log_accept_ratio + dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[1], log = TRUE) -
        dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
      if (i == 3) print('exp prior on')
    }

    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      gamma <- gamma_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }

    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma3 =  sigma3*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #*****************************************************
    #Beta-Gamma TRANSFORM
    if(FLAGS_LIST$ABG_TRANSFORM){
      gamma_dash <- gamma + rnorm(1, sd = sigma4)

      #Prior > 1 #* TRY WITHOUT REFLECTION
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash
      }
      #New b
       beta_transform = ((alpha + beta*gamma) - alpha)/gamma_dash #beta = (r0 - a)c

      if( beta_transform >= 0){ #Only accept values of beta> 0

        logl_new = LOG_LIKE_SSE_LSE(epidemic_data, alpha, beta_transform, gamma_dash)
        log_accept_ratio = logl_new - log_like

        #PRIORS
        #Beta prior
        if (FLAGS_LIST$BETA_PRIOR_GA) {
          tot_beta_prior = dgamma(beta_transform, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE) -
            dgamma(beta, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE)
        } else {
          tot_beta_prior = - beta_transform + beta #exp(1) prior
        }

        #gamma prior
        if (FLAGS_LIST$GAMMA_PRIOR_GA) {
          tot_gamma_prior = dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE) -
            dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
        } else {
          tot_gamma_prior = - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
        }

        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + tot_beta_prior + tot_gamma_prior

        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          beta <- beta_transform
          gamma <- gamma_dash
          log_like <- logl_new
          list_accept_counts$count_accept4 = list_accept_counts$count_accept4 + 1
        }

        #Sigma (Adpative)
        if (FLAGS_LIST$ADAPTIVE){
          accept_prob = min(1, exp(log_accept_ratio))
          sigma4 = sigma4*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    }

    #*****************************************************
    #Alpha-Gamma TRANSFORM
    if(FLAGS_LIST$ABG_TRANSFORM){

      gamma_dash <- gamma+ rnorm(1, sd = sigma5)
      #Prior > 1
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash
      }
      #New alpha
       alpha_transform = (alpha + beta*gamma) - beta*gamma_dash #alpha = (r0 - beta*gamma)

      if( alpha_transform >= 0){ #Only accept values of beta> 0

        logl_new = LOG_LIKE_SSE_LSE(epidemic_data, alpha_transform, beta, gamma_dash)
        log_accept_ratio = logl_new - log_like

        #PRIORS
        #gamma prior
        if (FLAGS_LIST$GAMMA_PRIOR_GA) {
          tot_gamma_prior = dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE) -
            dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
        } else {
          tot_gamma_prior = - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
        }

        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio - alpha_transform + alpha + tot_gamma_prior

        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          alpha <- alpha_transform
          gamma <- gamma_dash
          log_like <- logl_new
          list_accept_counts$count_accept5 = list_accept_counts$count_accept5 + 1
        }

        #Sigma (Adpative)
        if (FLAGS_LIST$ADAPTIVE){
          accept_prob = min(1, exp(log_accept_ratio))
          sigma5 = sigma5*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    }

    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      #print(paste0('i = ', i))
      i_thin = i/thinning_factor
      alpha_vec[i_thin] <- alpha; beta_vec[i_thin] <- beta
      gamma_vec[i_thin] <- gamma; r0_vec[i_thin] <- alpha + beta*gamma
      log_like_vec[i_thin] <- log_like
      sigma$sigma1_vec[i_thin] = sigma1; sigma$sigma2_vec[i_thin] = sigma2; sigma$sigma3_vec[i_thin] = sigma3
      sigma$sigma4_vec[i_thin] = sigma4; sigma$sigma5_vec[i_thin] = sigma5
    }
  }

  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1)
  accept_rate4 = 100*list_accept_counts$count_accept4/(n_mcmc-1)
  accept_rate5 = 100*list_accept_counts$count_accept5/(n_mcmc-1)

  #Acceptance rates
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5)
  print(list_accept_rates)

  #Return a, acceptance rate
  return(list(alpha_vec = alpha_vec, beta_vec = beta_vec, gamma_vec = gamma_vec, r0_vec = r0_vec,
              log_like_vec = log_like_vec, sigma = sigma,
              list_accept_rates = list_accept_rates))
}
