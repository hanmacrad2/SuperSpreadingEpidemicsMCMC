#********************************************************
#1. INDIVIDUAL R0 MCMC WITH ADAPTIVE SHAPING                           
#********************************************************

#*********************************************
#* SIMULATE ssi Model - Individual reproduction number
#**********************************************
#' @export
SIMULATE_EPI_SSI <- function(num_days = 50, R0 = 1.6, k = 1.1,
                        shape_gamma = 6, scale_gamma = 1,
                        epi_data = c(0,0,0), SIM_DATA = TRUE) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); 
  eta_vec = vector('numeric', num_days); 
  
  if(SIM_DATA){
    x[1] = 2
  } else {
    x[1] = epi_data[1] 
  }
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Infectiousness Pressure
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #ETA (t-1)
    eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = R0/k) #Draw eta from previous time step
    
    #INFECTIVITY
    infectivity = rev(prob_infect[1:t-1]) 
    
    #POISSON; OFFSPRING DISTRIBUTION
    #browser()
    total_rate = sum(eta_vec[1:t-1]*infectivity) 
    
    x[t] = rpois(1, total_rate)
    
  }
  return(list(epidemic_data = x, eta_vec = eta_vec))
}

#**********************************************
#LOG LIKELIHOOD
#**********************************************
#' @export
LOG_LIKE_SSI <- function(epi_data, infect_curve_ga, ssi_params, eta, FLAG_MCMC = TRUE){ #eta - a vector of length epi_data. eta[1] = infectivity of epi_datat[1]
  
  #Params
  num_days = length(epi_data)
  R0 = ssi_params[1]; k = ssi_params[2]
  loglike = 0;
  #eta = eta[eta > 0]
  
  for (t in 2:num_days) {
     
    #POISSON DENSITY
    total_poi_rate = sum(eta[1:(t-1)]*rev(infect_curve_ga[1:t-1]))
    
    if (total_poi_rate < 0) {
      log_poi_prob = -Inf
    } else {
      log_poi_prob = dpois(epi_data[t], lambda = total_poi_rate, log = TRUE) #NaN if eta == negative 
    }
    
    loglike = loglike + log_poi_prob 
    
    if (FLAG_MCMC && epi_data[t-1] > 0) {
      
      #ETA; GAMMA
      log_eta_prob = dgamma(eta[t-1], shape = epi_data[t-1]*k, scale = R0/k, log = TRUE)                  #dgamma with shape 0 == -Inf
      loglike = loglike + log_eta_prob
      
    } 
  }
  # 
  # if(loglike > 0 || is.nan(loglike)){
  #   browser()
  # }
  return(loglike)
}

#PRIOR
SET_SSI_PRIOR <- function(params, params_dash){
  
  #PRIORS
  list_priors = GET_LIST_PRIORS_SSI() 
  PRIORS_USED =  GET_PRIORS_USED() 
  
  R0 = params[1]; R0_dash = params_dash[1]
  k =  params[2]; k_dash = params_dash[2]
  
  if(PRIORS_USED$SSI$R0$EXP){
    
    p = dexp(R0_dash, rate = list_priors$r0[1], log = TRUE) -
      dexp(R0, rate = list_priors$r0[1], log = TRUE)
  }
  
  if (PRIORS_USED$SSI$k$EXP) {
    
    p = p + dexp(k_dash, rate = list_priors$k[1], log = TRUE) -
      dexp(k, rate = list_priors$k[1], log = TRUE)
  }
  
  return(p) 
}

    
#********************************************************
#1. MCMC INFERENCE FOR ssi MODEL - INDIVIDUAL R0  (INC. ADAPTIVE SCALING)                           
#********************************************************
#' @export
MCMC_INFER_SSI <- function(epidemic_data, n_mcmc = 50000,
                              mcmc_inputs = list(mod_start_points = c(1.2, 0.16),
                                                 dim = 2, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(R0_prior = c(1, 0), k_prior = c()),
                                                 thinning_factor = 10, burn_in_pc = 0.2),
                              FLAGS_LIST = list(BURN_IN = TRUE, ADAPTIVE = TRUE, THIN = TRUE)) {    
  
  #MCMC PARAMS + VECTORS
  num_days = length(epidemic_data); 
  vec_min = rep(0, mcmc_inputs$dim);
  count_accept = 0; count_accept_da = 0
  i_thin = 1
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #BURN_IN: Initial samples are not completely valid; the Markov Chain has not stabilized to the stationary distribution. The burn in samples allow you to discard these initial samples that are not yet at the stationary.
  if(FLAGS_LIST$BURN_IN){
    burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; print(paste0('N burn-in = ', burn_in_start))
    #Adjust mcmc vector store size
    mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size; ; print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
    mcmc_vec_size = mcmc_vec_size + 1
    }
  
  print(paste0('mcmc_vec_size', mcmc_vec_size))
  
  #MODEL PARAMETERS
  infect_curve_ga = GET_INFECT_GAMMA_CURVE(epidemic_data)
  #infectivity_vec = get_infectious_curve(epidemic_data)
  
  eta_dim = length(epidemic_data)-1 
  eta = epidemic_data[1:eta_dim];
  eta_matrix = matrix(NA, mcmc_vec_size, eta_dim); 
  
  ssi_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   #Changed from 0 to NA (As should be overwriting all cases)
  ssi_params_matrix[1,] <- mcmc_inputs$mod_start_points; ssi_params = ssi_params_matrix[1,] #2x1 #as.matrix
  
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSI(epidemic_data, infect_curve_ga, ssi_params, eta);  log_like = log_like_vec[1]
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  scaling_vec <- vector('numeric', mcmc_vec_size); scaling_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssi_params + ssi_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssi_params_matrix[1,]) + tcrossprod(ssi_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #SIGMA ETA (ONE DIMENSION - ROBINS MUNROE) 
  sigma_eta =  0.5*rep(1, eta_dim)
  sigma_eta_matrix = matrix(0, mcmc_vec_size, eta_dim); sigma_eta_matrix[1,] =  sigma_eta;

  #MCMC
  for(i in 2:n_mcmc) {

    if(i%%5000 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    ssi_params_dash = c(ssi_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssi_params_dash - vec_min) >= 0){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSI(epidemic_data, infect_curve_ga, ssi_params_dash, eta)
      log_accept_ratio = logl_new - log_like
      #PRIORS
      log_accept_ratio = log_accept_ratio + SET_SSI_PRIOR(ssi_params, ssi_params_dash)
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        ssi_params <- ssi_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssi_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssi_params)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
      
      #ACCEPTANCE PROB - LAMBDA ADAPTIVE SCALING
      accept_prob = min(1, exp(log_accept_ratio))
      
    } else {
      
      accept_prob = 0
    }
    
    #LAMBDA ADAPTIVE SCALING
    scaling =  scaling*exp(delta/i*(accept_prob - mcmc_inputs$target_acceptance_rate))
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    for(t in which(epidemic_data[1:(num_days-1)] > 0)){ #Only update eta's where x[t] > 0 epidemic_data[1:(num_days-1)]
      
      v = rep(0, length(eta))
      v[t] = 1
      
      #METROPOLIS STEP 
      eta_dash = abs(eta + rnorm(1, 0, sigma_eta[t])*v) #normalise the t_th element of eta #or variance = x[t]
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSI(epidemic_data, infect_curve_ga, ssi_params, eta_dash)
      log_accept_ratio = logl_new - log_like
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        
        #ACCEPT
        eta <- eta_dash
        log_like <- logl_new
        count_accept_da = count_accept_da + 1
      }
      
      accept_prob = min(1, exp(log_accept_ratio)) #Acceptance PROB = MIN(1, EXP(ACCPET_PROB))
      sigma_eta[t] =  sigma_eta[t]*exp(delta/(1+i)*(accept_prob - mcmc_inputs$target_acceptance_rate))
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
      
      ssi_params_matrix[i_thin,] = ssi_params
      log_like_vec[i_thin] <- log_like
      scaling_vec[i_thin] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      eta_matrix[i_thin, ] <- eta 
      sigma_eta_matrix[i_thin, ] = sigma_eta
      i_thin = i_thin + 1
    }
    
  } #END FOR LOOP

  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*num_days)
  
  #Return a, acceptance rate
  return(list(ssi_params_matrix = ssi_params_matrix, eta_matrix = eta_matrix,
              sigma_eta_matrix = sigma_eta_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate, accept_rate_da = accept_rate_da))
} 