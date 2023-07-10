#SSNB model
library(MASS)
#SIMULATE
#' @export
SIMULATE_EPI_SSNB <- function(num_days = 30, R0 = 1.6, k = 0.16,
                              shape_gamma = 6, scale_gamma = 1,
                              epi_data = c(0,0,0), SIM_DATA = TRUE,
                              FLAG_NEGBIN_PARAMATERISATION = list(param_mu = TRUE, param_prob = FALSE)) {
  
  'Simulate an epidemic with Superspreading events
  alpha - R0'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); 
  
  if(SIM_DATA){
    x[1] = 2
  } else {
    x[1] = epi_data[1] 
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    
    #NEGATIVE BINOMIAL PARAMETERISATION
    if (FLAG_NEGBIN_PARAMATERISATION$param_prob){
      x[t] = rnbinom(1, size = k*lambda_t, prob =  k/(R0 + k)) #Neg Bin parameterisation #1
    } else if (FLAG_NEGBIN_PARAMATERISATION$param_mu) {
      x[t] = rnbinom(1, size = k, mu =  R0*lambda_t) #Neg Bin parameterisation #2
    }
  }
  
  return(x)
}

#************************
#* LOG LIKELIHOOD SSNB
#* ***********************
#' @export
LOG_LIKE_SSNB <- function(epidemic_data, lambda_vec, ssnb_params){
  
  #Params
  k = ssnb_params[1]; R0 = ssnb_params[2] 
  num_days = length(epidemic_data); loglike = 0
  
  for (t in 2:num_days) {
    
    loglike_t = dnbinom(epidemic_data[t],
                        size = k, mu =  R0*lambda_vec[t], log = TRUE) #Neg Bin parameterisation #1 
    
    if(!is.nan(loglike_t)) { 
      
      loglike = loglike + loglike_t
      
    } else {
      print('log likelihood: NA')
      print(paste0('k, R0', k, R0))
    }
  }

  return(loglike)
  
}

#************************
#* LOG LIKELIHOOD SSNB
#* ***********************
#' @export
LOG_LIKE_SSNB_PARAMETERISATIONS <- function(x, lambda_vec, ssnb_params, 
                          FLAG_NEGBIN_PARAMATERISATION = list(param_mu = TRUE, param_prob = FALSE)){
  
  #Params
  k = ssnb_params[1]; R0 = ssnb_params[2]
  num_days = length(x); loglike = 0
  
  for (t in 2:num_days) {
    
    #NEGATIVE BINOMIAL PARAMETERISATION
    if (FLAG_NEGBIN_PARAMATERISATION$param_mu){
      
      loglike_t = dnbinom(x[t], size = k, mu =  R0*lambda_vec[t], log = TRUE) #Neg Bin parameterisation #1 
      
      if(!is.na(loglike_t)) { #likelihood = 0; log_likelihood = -Inf
        
        loglike = loglike + loglike_t 
      }
      
    } else if (FLAG_NEGBIN_PARAMATERISATION$param_prob) {
      loglike_t = dnbinom(x[t], size = k*lambda_vec[t], prob =  k/(R0 + k), log = TRUE) #Neg Bin parameterisation #2
      
      if(!is.na(loglike_t)) { 
        loglike = loglike + loglike_t 
      }
    }
  }
  
  return(loglike)
}

#********************************************************
#1. MCMC INFERENCE FOR SSIC MODEL - INDIVIDUAL R0  (INC. ADAPTIVE SCALING)                           
#********************************************************
#' @export
MCMC_INFER_SSNB <- function(epidemic_data, n_mcmc,
                            mcmc_inputs = list(mod_start_points = c(0.16, 1.8),
                                               dim = 2, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                               thinning_factor = 10, burn_in_pc = 0.2),
                            priors = list(pk_exp = c(1,0), pk_ga_shape = 0.001, pk_ga_rte = 0.001, 
                                          pr0_unif = c(1.0,4), p_prob_unif = c(0,1)),
                            PRIORS_USED = list(EXP_K = TRUE, GAMMA_K = FALSE),
                            FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE, PRIOR = TRUE, BURN_IN = TRUE),
                            FLAG_NEGBIN_PARAMATERISATION = list(param_mu = TRUE, param_prob = FALSE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper); #NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  print(FLAG_NEGBIN_PARAMATERISATION)
  num_days = length(epidemic_data)
  vec_min = rep(0, mcmc_inputs$dim)
  count_accept = 0; idx_thinned = 0
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #BURN_IN: 
  if(FLAGS_LIST$BURN_IN){
    burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; print(paste0('N burn-in = ', burn_in_start))
    #Adjust mcmc vector store size
    mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size; ; print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
    #mcmc_vec_size = ceil(mcmc_vec_size)
  }
  
  #MODEL PARAMETERS
  lambda_vec = get_lambda(epidemic_data)
  ssnb_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   #Changed from 0 to NA (As should be overwriting all cases)
  ssnb_params_matrix[1,] <- mcmc_inputs$mod_start_points; ssnb_params = ssnb_params_matrix[1,] #2x1 #as.matrix
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSNB(epidemic_data, lambda_vec, ssnb_params) #, FLAG_NEGBIN_PARAMATERISATION)
  log_like = log_like_vec[1]
  #pk_ga_scale = ((priors$negbin_k_prior_ga_sd)^2)/priors$negbin_k_prior_ga_mean
  #pk_ga_shape = negbin_scale*priors$negbin_k_prior_ga_mean
  print(paste0('PRIOR ON K', PRIORS_USED))
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  scaling_vec <- vector('numeric', mcmc_vec_size); scaling_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssnb_params + ssnb_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssnb_params_matrix[1,]) + tcrossprod(ssnb_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #MCMC
  for(i in 2:n_mcmc) {

    if(i%%10000 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    #print(paste0('scaling*c_star*sigma_i', scaling*c_star*sigma_i))
    ssnb_params_dash = c(ssnb_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssnb_params_dash - vec_min) >= 0){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSNB(epidemic_data, lambda_vec, ssnb_params_dash) #, FLAG_NEGBIN_PARAMATERISATION)
      
      #ACCEPTANCE RATIO
      log_accept_ratio = logl_new - log_like
      
      #PRIORS
      #EXTRACT PARAMS FORPRIORS
      k =  ssnb_params[1]; k_dash = ssnb_params_dash[1]
      R0 = ssnb_params[2]; R0_dash = ssnb_params_dash[2]
    
      if(FLAG_NEGBIN_PARAMATERISATION$param_mu && PRIORS_USED$EXP_K) {

        log_accept_ratio = log_accept_ratio +
          dexp(k_dash, rate = priors$pk_exp[1], log = TRUE) -
          dexp(k, rate = priors$pk_exp[1], log = TRUE) +
        dunif(R0_dash, min = priors$pr0_unif[1], max = priors$pr0_unif[2], log = TRUE) - #Uniform prior
        dunif(R0, min = priors$pr0_unif[1], max = priors$pr0_unif[2], log = TRUE)
        #priors_list$pr0_unif[1]*ssnb_params_dash[2] + priors_list$pr0_unif[1]*ssnb_params[2]

      } else if (FLAG_NEGBIN_PARAMATERISATION$param_mu && PRIORS_USED$GAMMA_K) {
        
        log_accept_ratio = log_accept_ratio +
          dgamma(k_dash, shape =  priors$pk_ga_shape, rate = priors$pk_ga_rte, log = TRUE) -
          dgamma(k, shape =  priors$pk_ga_shape,  rate = priors$pk_ga_rte, log = TRUE) +
        dunif(R0_dash, min = priors$pr0_unif[1], max = priors$pr0_unif[2], log = TRUE) - #Uniform prior
          dunif(R0, min = priors$pr0_unif[1], max = priors$pr0_unif[2], log = TRUE)
        
      } else if (FLAG_NEGBIN_PARAMATERISATION$param_prob && PRIORS_USED$GAMMA_K){

        log_accept_ratio = log_accept_ratio +
          dgamma(k_dash, shape =  priors$pk_ga_shape, rate = priors$pk_ga_rte, log = TRUE) -
          dgamma(k, shape =  priors$pk_ga_shape,  rate = priors$pk_ga_rte, log = TRUE) +
          dunif(k_dash/(R0_dash + k_dash), log = TRUE) -  dunif(k/(R0 + k), log = TRUE) #+
        #dunif(R0_dash, min = 0, max = 10 log = TRUE) -  dunif(R0, min = 0, max = 10 log = TRUE)
      }
 
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        ssnb_params <- ssnb_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssnb_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssnb_params)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
      
      #ACCEPTANCE PROB - ADAPTIVE SCALING
      accept_prob = min(1, exp(log_accept_ratio))
      
    } else {
     accept_prob = 0
    }
    
    #ADAPTIVE SCALING (needs acceptance probability)
    scaling =  scaling*exp(delta/i*(accept_prob - mcmc_inputs$target_acceptance_rate))
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i > burn_in_start & i%%thinning_factor == 0) {
      #print(paste0('i = ', i))
      #i_thin = i/thinning_factor; 
      #print(paste0('i = ', i))
      idx_thinned = idx_thinned + 1
      #print(paste0('idx_thinned = ', idx_thinned))
      ssnb_params_matrix[idx_thinned,] = ssnb_params
      log_like_vec[idx_thinned] <- log_like
      scaling_vec[idx_thinned] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
    }
    
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  
  #Return a, acceptance rate
  return(list(ssnb_params_matrix = ssnb_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate))
} 