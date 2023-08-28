#SSE model

#SIMULATE
#' @export
SIMULATE_EPI_SSE <- function(num_days = 30, R0 = 1.6, k = 0.16,
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
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    
    #NEGATIVE BINOMIAL PARAMETERISATION
    if (FLAG_NEGBIN_PARAMATERISATION$param_prob){
      
      x[t] = rnbinom(1, size = k*lambda_t, prob =  k/(R0 + k)) #Neg Bin parameterisation #1
    
      } else if (FLAG_NEGBIN_PARAMATERISATION$param_mu) {
      
      x[t] = rnbinom(1, size = k*lambda_t, mu =  R0*lambda_t)
      #x[t] = rnbinom(1, size = k, mu =  R0*lambda_t) #Neg Bin parameterisation #2
    }
  }
  
  return(x)
}

#************************
#* LOG LIKELIHOOD SSE
#* ***********************
#' @export
LOG_LIKE_SSE <- function(epidemic_data, lambda_vec, sse_params){
  
  #Params
  R0 = sse_params[1];  k = sse_params[2]
  num_days = length(epidemic_data); loglike = 0
  
  for (t in 2:num_days) {
    
    loglike_t = dnbinom(epidemic_data[t],
                        size = k*lambda_vec[t], mu =  R0*lambda_vec[t], log = TRUE) #Neg Bin parameterisation #1 
    
    if(!is.nan(loglike_t)) { 
      
      loglike = loglike + loglike_t
      
    } else {
      print('log likelihood: NaN')
      print(paste0('k, R0', k, R0))
    }
  }

  return(loglike)
  
}

#************************
#* PRIOR SSE
#* ***********************
SET_SSE_PRIOR <- function(params, params_dash){
  
  #PRIORS
  list_priors = GET_LIST_PRIORS_SSE() 
  PRIORS_USED =  GET_PRIORS_USED() 
  
  R0 = params[1]; R0_dash = params_dash[1]
  k =  params[2]; k_dash = params_dash[2]
  
  if(PRIORS_USED$SSE$R0$EXP && PRIORS_USED$SSE$k$EXP){
    
      p = dexp(R0_dash, rate = list_priors$r0[1], log = TRUE) -
        dexp(R0, rate = list_priors$r0[1], log = TRUE) +
        dexp(k_dash, rate = list_priors$k[1], log = TRUE) -
        dexp(k, rate = list_priors$k[1], log = TRUE)
  }
 
  return(p) 
}


#********************************************************
#1. MCMC INFERENCE FOR SSIC MODEL - INDIVIDUAL R0  (INC. ADAPTIVE SCALING)                           
#********************************************************
#' @export
MCMC_INFER_SSE <- function(epidemic_data, n_mcmc,
                            mcmc_inputs = list(mod_start_points = c(1.2, 0.5),
                                               dim = 2, target_acceptance_rate = 0.4, v0 = 100,  
                                               thinning_factor = 10, burn_in_pc = 0.2),
                            FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE, BURN_IN = TRUE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper); #NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  num_days = length(epidemic_data)
  vec_min = rep(0, mcmc_inputs$dim)
  count_accept = 0; i_thin = 1
  
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
  sse_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);  
  sse_params_matrix[1,] <- mcmc_inputs$mod_start_points; sse_params = sse_params_matrix[1,] 
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSE(epidemic_data, lambda_vec, sse_params) 
  log_like = log_like_vec[1]

  
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

      #Priors
      log_accept_ratio = log_accept_ratio + SET_SSE_PRIOR(sse_params, sse_params_dash)

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
      scaling_vec[i_thin] <- scaling 
      i_thin = i_thin + 1                 #Taking role of sigma, overall scaling constant
                                          #Sigma becomes estimate of the covariance matrix of the posterior
    }
    
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  
  #Return a, acceptance rate
  return(list(sse_params_matrix = sse_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate))
} 
