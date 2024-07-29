#SSE model

#SIMULATE
#' @export
SIMULATE_EPI_SSE <- function(num_days = 50, r0 = 2.0, k = 0.2, #k = 0.16, 0.8
                             shape_gamma = 6, scale_gamma = 1,
                             epi_data = c(0,0,0), SIM_DATA = TRUE,
                             FLAG_NEGBIN_PARAMATERISATION = list(param_mu = TRUE, param_prob = FALSE)) {
  
  'Simulate an epidemic with Superspreading events
  alpha - r0'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); 
  
  if(SIM_DATA){
    x[1] = 2
  } else {
    x[1] = epi_data[1] 
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    
    #NEGATIVE BINOMIAL PARAMETERISATION
    if (FLAG_NEGBIN_PARAMATERISATION$param_mu){
      
      x[t] = rnbinom(1, size = k*lambda_t, mu =  r0*lambda_t) #Neg Bin parameterisation #1
      
    } else if (FLAG_NEGBIN_PARAMATERISATION$param_prob) {
      
      x[t] = rnbinom(1, size = k*lambda_t, prob =  k/(r0 + k)) #Neg Bin parameterisation #2
      
    }
  }
  
  return(x)
}

#************************
#* LOG LIKELIHOOD sse
#* ***********************
#' @export
LOG_LIKE_SSE <- function(epidemic_data, lambda_vec, sse_params){
  
  #Params
  r0 = sse_params[1];  k = sse_params[2]
  num_days = length(epidemic_data); loglike = 0
  
  for (t in 2:num_days) {
    
    loglike_t = dnbinom(epidemic_data[t],
                        size = k*lambda_vec[t], mu =  r0*lambda_vec[t], log = TRUE) #Neg Bin parameterisation #1 
    
    if(!is.nan(loglike_t)) { 
      
      loglike = loglike + loglike_t
      
    } else {
      print('log likelihood: NaN')
      print(paste0('k, r0', k, r0))
    }
  }
  
  return(loglike)
  
}

#************************
#* LOG LIKELIHOOD sse
#* ***********************
#' @export
LOG_LIKE_SSE_POINTWISE <- function(epidemic_data, lambda_vec, sse_params){
  
  #Params
  r0 = sse_params[1];  k = sse_params[2]
  num_days = length(epidemic_data); loglike = 0
  log_like_vec = length(epidemic_data-1)
  
  for (t in 2:num_days) {
    
    log_like_vec[t-1] = dnbinom(epidemic_data[t],
                        size = k*lambda_vec[t], mu =  r0*lambda_vec[t], log = TRUE) #Neg Bin parameterisation #1 
    
    
  }
  
  return(log_like_vec)
  
}

#****************
#* SET PRIOR SSE
#* **************
SET_SSE_PRIOR <- function(sse_params, sse_params_dash, PRIORS_USED){
  
  #EXTRACT PARAMS FOR PRIORS
  r0 = sse_params[1]; r0_dash = sse_params_dash[1]
  k =  sse_params[2]; k_dash = sse_params_dash[2]
  list_priors = GET_LIST_PRIORS_SSE() 
  prior = 0
  
  if (PRIORS_USED$SSE$r0$EXP) {
    
    prior = dexp(r0_dash, rate = list_priors$r0[1], log = TRUE) -
      dexp(r0, rate = list_priors$r0[1], log = TRUE) 
    
  } else if (PRIORS_USED$SSE$r0$GAMMA){
    #print('Gamma prior used r0 sse')
    
    prior = dgamma(r0_dash, shape = list_priors$r0_gamma[1], scale = list_priors$r0_gamma[2], log = TRUE) -
      dgamma(r0, shape = list_priors$r0_gamma[1], scale = list_priors$r0_gamma[2], log = TRUE) 
  
    } else if (PRIORS_USED$SSE$r0$UNIF){
    
      prior = dunif(r0_dash, min = list_priors$r0_unif[1], max = list_priors$r0_unif[2], log = TRUE) -
        dunif(r0, min = list_priors$r0_unif[1], max = list_priors$r0_unif[2], log = TRUE) 
  }
  
  if (PRIORS_USED$SSE$k$EXP) {
    prior = prior + dexp(k_dash, rate = list_priors$k[1], log = TRUE) -
      dexp(k, rate = list_priors$k[1], log = TRUE) 
    
  } else if (PRIORS_USED$SSE$k$UNIF) {
    
    prior = prior + dunif(k_dash, min = list_priors$k_unif[1], max = list_priors$k_unif[2], log = TRUE) -
      dunif(k, min = list_priors$k_unif[1], max = list_priors$k_unif[2], log = TRUE)
  }
  
  
  return(prior)
}

#********************************************************
#1. MCMC INFERENCE FOR SSE MODEL  #24/12/23 param_starts = c(1.2, 0.15), #= GET_PRIORS_USED(), 
#********************************************************
#' @export
MCMC_INFER_SSE <- function(epidemic_data, n_mcmc, PRIORS_USED = GET_PRIORS_USED(),
                           param_starts = c(2.0, 0.2), #Changed 08/01/24 c(1.0, 0.5)
                           mcmc_inputs = list(dim = 2, target_acceptance_rate = 0.234, #0.4, 
                                              v0 = 100, thinning_factor = 10, burn_in_pc = 0.2),
                           FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE,
                                             PRIOR = TRUE, BURN_IN = TRUE, COMPUTE_WAIC = TRUE)){    
  
  #MCMC INITIAL POINTS
  r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  param_starts[1] = r0_start
  
  #MCMC PARAMS + VECTORS
  #print(paste0('PRIORS USED;', PRIORS_USED))
  print(paste0('PARAM STARTING POINTS;', param_starts))
  i_thin = 1
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
    mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size; print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
  }
  
  #MODEL PARAMETERS
  lambda_vec = get_lambda(epidemic_data)
  sse_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   #Changed from 0 to NA (As should be overwriting all cases)
  sse_params_matrix[1,] <- param_starts; sse_params = sse_params_matrix[1,] #2x1 #as.matrix
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSE(epidemic_data, lambda_vec, sse_params) #, FLAG_NEGBIN_PARAMATERISATION)
  log_like = log_like_vec[1]
  loglike_pointwise_matrix = matrix(nrow = mcmc_vec_size, ncol = length(epidemic_data)-1)
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  scaling_vec <- vector('numeric', mcmc_vec_size); scaling_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(sse_params + sse_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(sse_params_matrix[1,]) + tcrossprod(sse_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #MCMC
  for(i in 2:n_mcmc) {
    
    if(i%%10000 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    sse_params_dash = c(sse_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(sse_params_dash - vec_min) >= 0){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSE(epidemic_data, lambda_vec, sse_params_dash) #, FLAG_NEGBIN_PARAMATERISATION)
      
      #ACCEPTANCE RATIO
      log_accept_ratio = logl_new - log_like
      
      #PRIORS
      log_accept_ratio = log_accept_ratio + SET_SSE_PRIOR(sse_params, sse_params_dash, PRIORS_USED)
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        sse_params <- sse_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*sse_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(sse_params)
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
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
      sse_params_matrix[i_thin,] = sse_params
      log_like_vec[i_thin] <- log_like
      scaling_vec[i_thin] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      
      if (FLAGS_LIST$COMPUTE_WAIC){
        ll_vec = LOG_LIKE_SSE_POINTWISE(epidemic_data, lambda_vec, sse_params)
        loglike_pointwise_matrix[i_thin,] = ll_vec
      }
      
      i_thin = i_thin + 1
      
    }
    
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  
  #COMPUTE WAIC, DIC
  waic_result = WAIC(loglike_pointwise_matrix)$WAIC
  dic_result = GET_DIC(loglike_pointwise_matrix)
  
  #Return a, acceptance rate
  return(list(sse_params_matrix = sse_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              waic_result = waic_result, dic_result = dic_result,
              accept_rate = accept_rate, r0_start = r0_start))
} 


