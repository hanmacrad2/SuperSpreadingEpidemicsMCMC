#' Log likelihood of the Super-Spreading Individuals (SSI) epidemic model
#'
#' Returns the Log likelihood of the super-spreading individuals (SSI)epidemic model for given epidemic \code{"data"} and set of model parameters \code{"a, b"} and \code{"c"} and
#'
#' @param data data from the epidemic, namely daily infection counts
#' @param aX SSI model parameter \code{"a"}
#' @param bX SSI model parameter \code{"b"}
#' @param cX SSI model parameter \code{"c"}
#' @return Log likelihood of the SSI model
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' log_like = LOG_LIKE_SSI(epidemic_data, 0.8, 0.02, 20)
#'
LOG_LIKE_SSI <- function(sim_data, aX, bX, cX){

  #Data
  n = sim_data[[1]]; s = sim_data[[2]]

  #Params
  num_days = length(n)
  shape_gamma = 6; scale_gamma = 1
  logl = 0

  #INFECTIOUSNESS  - Difference of 2 GAMMA distributions. Discretized
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)

  for (t in 1:num_days) { #*1 or 2

    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS
    lambda_t = sum((n[1:(t-1)] + cX*s[1:(t-1)])*rev(prob_infect[1:(t-1)]))

    #LOG-LIKELIHOOD
    logl = logl - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
  }

  logl
}

#' MCMC adaptive algorithm for Super-Spreading Individuals epidemic model
#'
#' MCMC algorithm with Adaptation for obtaining samples from the parameters of a
#' super-spreading individuals (SSI) epidemic model. Includes; Adaptation, Data Augmentation, B-C & A-C transform
#'
#' @param data data from the epidemic, namely daily infection counts
#' @param mcmc_inputs A list of mcmc specifications including
#' \itemize{
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
#'   \item \code{"mod_start_points"} - Model parameter starting points; where the mcmc algorithm begins sampling from
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
#'mcmc_inputs = list(n_mcmc = 500000, mod_start_points = list(m1 = 0.72, m2 = 0.0038, m3 = 22),
#' alpha_star = 0.4, thinning_factor = 10)
#'
#' #START MCMC
#'mcmc_ssi_output = SSI_MCMC_ADAPTIVE(epidemic_data, mcmc_inputs)
#'

#************************************************************************
#1. SSI MCMC                              (W/ DATA AUGMENTATION OPTION)
#************************************************************************
SSI_MCMC_ADAPTIVE <- function(data,
                              mcmc_inputs = list(n_mcmc = 500000,
                                                 mod_start_points = list(m1 = 0.72, m2 = 0.0038, m3 = 22), alpha_star = 0.4,
                                                 thinning_factor = 10),
                              priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0),
                                                 c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                              FLAGS_LIST = list(ADAPTIVE = TRUE, DATA_AUG = TRUE, BCA_TRANSFORM = TRUE,
                                                PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE,
                                                THIN = TRUE)) {

  'Returns MCMC samples of SSI model parameters (a, b, c, r0 = a + b*c)
  w/ acceptance rates.
  INCLUDES; ADAPTATION, DATA AUGMENTATION, B-C & A-C transform'

  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'


  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  time = length(data[[1]]); n_mcmc = mcmc_inputs$n_mcmc;
  print(paste0('num mcmc iters = ', n_mcmc))

  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1
    mcmc_vec_size = n_mcmc
  }

  #INITIALISE MCMC VECTORS
  a_vec <- vector('numeric', mcmc_vec_size); b_vec <- vector('numeric', mcmc_vec_size)
  c_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);

  #INITIALISE MCMC[1]
  a_vec[1] <- mcmc_inputs$mod_start_points$m1; b_vec[1] <- mcmc_inputs$mod_start_points$m2
  c_vec[1] <- mcmc_inputs$mod_start_points$m3; r0_vec[1] <- a_vec[1] + b_vec[1]*c_vec[1]
  log_like_vec[1] <- LOG_LIKE_SSI(data, a_vec[1], b_vec[1], c_vec[1])

  #INITIALISE RUNNING PARAMS
  a = a_vec[1]; b = b_vec[1]; c = c_vec[1]; log_like = log_like_vec[1]

  #SIGMA
  sigma1 =  0.4*mcmc_inputs$mod_start_points$m1;  sigma2 = 0.3*mcmc_inputs$mod_start_points$m2
  sigma3 = 0.5*mcmc_inputs$mod_start_points$m3; sigma4 = 0.85*mcmc_inputs$mod_start_points$m3
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
    sigma = list(sigma1 <- sigma1, sigma2 <- sigma2,
                 sigma3 <- sigma3, sigma4 <- sigma4,
                 sigma5 <- sigma5)
  }

  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0, count_accept6 = 0)

  #DATA AUG OUTPUT
  #mat_count_da = matrix(0, mcmc_vec_size, time) #i x t
  non_ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR

  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {

    #****************************************************** s
    #a
    a_dash <- a + rnorm(1, sd = sigma1)

    if(a_dash < 0){
      a_dash = abs(a_dash)
    }

    #log a
    logl_new = LOG_LIKE_SSI(data, a_dash, b, c)
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_ratio = log_accept_ratio - a_dash + a #*Actually this is the Acceptance RATIO. ACCEPTANCE PROB = MIN(1, EXP(ACCPET_PROB))
    }

    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      a <- a_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    }

    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma1 =  sigma1*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #************************************************************************ Only if (b > 0){ ?
    #b
    b_dash <- b + rnorm(1, sd = sigma2)
    if(b_dash < 0){
      b_dash = abs(b_dash)
    }

    #loglikelihood
    logl_new = LOG_LIKE_SSI(data, a, b_dash, c)
    log_accept_ratio = logl_new - log_like

    #Priors
    if (FLAGS_LIST$B_PRIOR_GAMMA){
      log_accept_ratio = log_accept_ratio +
        dgamma(b_dash, shape = priors_list$b_prior_ga[1], scale = priors_list$b_prior_ga[2], log = TRUE) -
        dgamma(b, shape = priors_list$b_prior_ga[1], scale = priors_list$b_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - b_dash + b
    }

    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      b <- b_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    }

    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma2 =  sigma2*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #************************************************************************
    #c
    c_dash <- c + rnorm(1, sd = sigma3)
    if(c_dash < 1){
      c_dash = 2 - c_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSI(data, a, b, c_dash)
    log_accept_ratio = logl_new - log_like

    #Priors
    if(FLAGS_LIST$C_PRIOR_GAMMA){
      log_accept_ratio = log_accept_ratio + dgamma(c_dash, shape = priors_list$c_prior_ga[1], scale = priors_list$c_prior_ga[1], log = TRUE) -
        dgamma(c, shape = priors_list$c_prior_ga[1], scale = priors_list$c_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c
      if (i == 3) print('exp prior on')
    }

    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      c <- c_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }

    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma3 =  sigma3*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #*****************************************************
    #B-C TRANSFORM
    if(FLAGS_LIST$BCA_TRANSFORM){

      c_dash <- c + rnorm(1, sd = sigma4)
      #Prior > 1 #* TRY WITHOUT REFLECTION
      if(c_dash < 1){
        c_dash = 2 - c_dash
      }
      #New b
      b_transform = ((a + b*c) - a)/c_dash #b = (r0 - a)c

      if(b_transform >= 0){ #Only accept values of b > 0

        logl_new = LOG_LIKE_SSI(data, a, b_transform, c_dash)
        log_accept_ratio = logl_new - log_like

        #PRIORS
        #b prior
        if (FLAGS_LIST$B_PRIOR_GAMMA) {
          tot_b_prior = dgamma(b_transform, shape = priors_list$b_prior_ga[1], scale = priors_list$b_prior_ga[2], log = TRUE) -
            dgamma(b, shape = priors_list$b_prior_ga[1], scale = priors_list$b_prior_ga[2], log = TRUE)
        } else {
          tot_b_prior = - b_transform + b #exp(1) piror
        }

        #c prior
        if (FLAGS_LIST$C_PRIOR_GAMMA) {
          tot_c_prior = dgamma(c_dash, shape = priors_list$c_prior_ga[1], scale = priors_list$c_prior_ga[2], log = TRUE) -
            dgamma(c, shape = priors_list$c_prior_ga[1], scale = priors_list$c_prior_ga[2], log = TRUE)
        } else {
          tot_c_prior = - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c
        }

        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + tot_b_prior + tot_c_prior

        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          b <- b_transform
          c <- c_dash
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
    #A-C TRANSFORM
    if(FLAGS_LIST$BCA_TRANSFORM){

      c_dash <- c + rnorm(1, sd = sigma5)
      #Prior > 1
      if(c_dash < 1){
        c_dash = 2 - c_dash
      }
      #New a
      a_transform = (a + b*c) - b*c_dash #a = (r0 - b*c

      if(a_transform >= 0){ #Only accept values of b > 0

        logl_new = LOG_LIKE_SSI(data, a_transform, b, c_dash)
        log_accept_ratio = logl_new - log_like

        #PRIORS
        #c prior
        if (FLAGS_LIST$C_PRIOR_GAMMA) {
          tot_c_prior = dgamma(c_dash, shape = priors_list$c_prior_ga[1], scale = priors_list$c_prior_ga[2], log = TRUE) -
            dgamma(c, shape = priors_list$c_prior_ga[1], scale = priors_list$c_prior_ga[2], log = TRUE)
        } else {
          tot_c_prior = - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c
        }

        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + - a_transform + a + tot_c_prior

        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          a <- a_transform
          c <- c_dash
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

    #************************************
    #DATA AUGMENTATION
    #************************************
    if (FLAGS_LIST$DATA_AUG){

      #FOR EACH S_T
      for(t in 1:time){

        #Copy of data (or update as necessary)
        data_dash = data

        #STOCHASTIC PROPOSAL for s
        if (runif(1) < 0.5) {
          st_dash = data[[2]][t] + 1
        } else {
          st_dash = data[[2]][t] - 1
        }

        #CHECK
        if (st_dash < 0) {
          st_dash = 0
        }

        #ACCEPTANCE PROBABILITY
        data_dash[[2]][t] = st_dash #s_t = st_dash
        data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t

        if (data_dash[[1]][t] < 0){
          data_dash[[1]][t] = 0
        }
        #CRITERIA FOR S_T & N_T
        if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
          #Store
          non_ss[i/thinning_factor, t] = data[[1]][t]
          ss[i/thinning_factor, t] = data[[2]][t]
          next
        }

        logl_new = LOG_LIKE_SSI(data_dash, a, b, c)
        log_accept_ratio = logl_new - log_like

        #METROPOLIS ACCEPTANCE STEP
        if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {

          #ACCEPT
          data <- data_dash
          log_like <- logl_new
          #mat_count_da[i, t] = mat_count_da[i, t] + 1
          list_accept_counts$count_accept6 = list_accept_counts$count_accept6 + 1
        }

        #Store
        non_ss[i/thinning_factor, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i/thinning_factor, t] = data[[2]][t]
      }
    }

    #Loglikelihood Check (Passing - no error)
    #if (!(log_like && LOG_LIKE_SSI(data, a, b, c))) print('loglike doesnt exist')
    #else (log_like!=LOG_LIKE_SSI(data, a, b, c)) print(paste0('ERROR! logl diff = ', log_like - LOG_LIKE_SSI(data, a, b, c)))

    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      #print(paste0('i = ', i))
      i_thin = i/thinning_factor
      a_vec[i_thin] <- a; b_vec[i_thin] <- b
      c_vec[i_thin] <- c; r0_vec[i_thin] <- a + b*c
      log_like_vec[i_thin] <- log_like #PLOT!!
      sigma$sigma1[i_thin] = sigma1; sigma$sigma2[i_thin] = sigma2; sigma$sigma3[i_thin] = sigma3
      sigma$sigma4[i_thin] = sigma4; sigma$sigma5[i_thin] = sigma5
      }
  }

  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1)
  accept_rate4 = 100*list_accept_counts$count_accept4/(n_mcmc-1)
  accept_rate5 = 100*list_accept_counts$count_accept5/(n_mcmc-1)
  accept_rate6 = 100*list_accept_counts$count_accept6/((n_mcmc-1)*time) #i x t

  #Acceptance rates
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5,
                           accept_rate6 = accept_rate6)
  print(list_accept_rates)

  #Return a, acceptance rate
  return(list(a_vec = a_vec, b_vec = b_vec, c_vec = c_vec, r0_vec = r0_vec,
              log_like_vec = log_like_vec, sigma = sigma,
              list_accept_rates = list_accept_rates,
              data = data, non_ss = non_ss, ss = ss))
}

