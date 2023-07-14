#********************************************************
#1. INDIVIDUAL R0 MCMC WITH ADAPTIVE SHAPING                           
#********************************************************
library(MASS)

#*********************************************
#* SIMULATE SSIR Model - Individual reproduction number
#**********************************************
#' @export
SIMULATE_EPI_SSIR <- function(num_days = 30, R0X = 1.6, k = 0.16,
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
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #ETA (t-1)
    eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = R0X/k) #Draw eta from previous time step
    
    #INFECTIVITY
    infectivity = rev(prob_infect[1:t]) 
    
    #POISSON; OFFSPRING DISTRIBUTION
    total_rate = sum(eta_vec*infectivity) #DOT PRODUCT
    
    x[t] = rpois(1, total_rate)
    
  }
  return(list(epidemic_data = x, eta_vec = eta_vec))
}

#**********************************************
#LOG LIKELIHOOD
#**********************************************
#' @export
LOG_LIKE_SSIR <- function(epi_data, infect_curve_ga, ssir_params, eta){ #eta - a vector of length epi_data. eta[1] = infectivity of epi_datat[1]
  
  #Params
  num_days = length(epi_data)
  R0 = ssir_params[1]; k = ssir_params[2]
  loglike = 0; count_na = 0; count_not_na = 0
  
  for (t in 2:num_days) {
    
    #POISSON
    total_poi_rate = sum(eta[1:(t-1)]*rev(infect_curve_ga[1:t-1]))
    log_poi_prob = dpois(epi_data[t], lambda = total_poi_rate, log = TRUE)
    
    if (!is.nan(log_poi_prob) && !is.na(log_poi_prob) && !is.infinite(log_poi_prob)){
      loglike = loglike + log_poi_prob 
      count_not_na = count_not_na +1
      #print(paste0('loglike: ', loglike))
      
    } else {
      count_na = count_na + 1
      #print('log_poi_prob == NaN')
      # print(paste0('log_poi_prob: ', log_poi_prob))
      # print(paste0('epi_data[t]: ', epi_data[t]))
      # print(paste0('eta[t-1]: ', eta[t-1]))
      # print(paste0('R0: ', R0))
      # print(paste0('k: ', k))
    }
    
    if (epi_data[t-1] > 0) {
      
      # print(paste0('t', t))
      # print(paste0('epi_data[t-1]: ', epi_data[t-1]))
      # print(paste0('eta[t-1]: ', eta[t-1]))
      # print(paste0('eta: ', eta))
      
      #ETA; GAMMA
      log_eta_prob = dgamma(eta[t-1], shape = epi_data[t-1]*k, scale = R0/k, log = TRUE) #dgamma with shape 0 == -Inf
     
       if (!is.nan(log_eta_prob) && !is.na(log_eta_prob) && !is.infinite(log_eta_prob)){
        
        loglike = loglike + log_eta_prob 
        count_not_na = count_not_na +1
        # print(paste0('t: ', t))
        # print(paste0('loglike: ', loglike))
      }
      
      else {
        count_na = count_na + 1
        # print(paste0('log_eta_prob: ', log_eta_prob))
        # # print('log_eta_prob == NaN')
        # print(paste0('epi_data[t]: ', epi_data[t]))
        # print(paste0('eta[t-1]: ', eta[t-1]))
        # print(paste0('R0: ', R0))
        # print(paste0('k: ', k))
      }
      
    } 
    
    # if (epi_data[t] == 0){
    #   print(paste0('epi_data[t]: ', epi_data[t]))
    #   print(paste0('log_poi_prob: ', log_poi_prob))
    #   print(paste0('log_eta_prob: ', log_eta_prob))
    # }
  }
  
  # print(paste0('loglike: ', loglike))
  # print(paste0('count_not_na', count_not_na))
  # print(paste0('count na', count_na))
  return(loglike)
}

#********************************************************
#1. MCMC INFERENCE FOR ssir MODEL - INDIVIDUAL R0  (INC. ADAPTIVE SCALING)                           
#********************************************************
#' @export
MCMC_INFER_SSIR <- function(epidemic_data, n_mcmc = 50000,
                              mcmc_inputs = list(mod_start_points = c(1.2, 0.16),
                                                 dim = 2, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(R0_prior = c(1, 0), k_prior = c()),
                                                 thinning_factor = 10, burn_in_pc = 0.2),
                              priors_list = list(R0_prior = c(1,0), k_prior = c(1, 0)),
                              FLAGS_LIST = list(BURN_IN = TRUE, ADAPTIVE = TRUE, THIN = TRUE, PRIOR_K1 = TRUE,
                                                PRIOR_K2 = FALSE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper); #NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
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
  
  eta_dim = length(epidemic_data)-1 #epidemic_data[1:(length(epidemic_data)-1)]
  
  eta = epidemic_data[1:eta_dim]; #print(paste0('eta = ', eta))
  eta_matrix = matrix(NA, mcmc_vec_size, eta_dim); 
  
  ssir_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   #Changed from 0 to NA (As should be overwriting all cases)
  ssir_params_matrix[1,] <- mcmc_inputs$mod_start_points; ssir_params = ssir_params_matrix[1,] #2x1 #as.matrix
  
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSIR(epidemic_data, infect_curve_ga, ssir_params, eta);  log_like = log_like_vec[1]
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  scaling_vec <- vector('numeric', mcmc_vec_size); scaling_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssir_params + ssir_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssir_params_matrix[1,]) + tcrossprod(ssir_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #SIGMA ETA (ONE DIMENSION - ROBINS MUNROE) 
  sigma_eta =  0.5*rep(1, eta_dim)
  sigma_eta_matrix = matrix(0, mcmc_vec_size, eta_dim); sigma_eta_matrix[1,] =  sigma_eta;
  
  #DIRECTORY - SAVING
  ifelse(!dir.exists(file.path(OUTER_FOLDER)), dir.create(file.path(OUTER_FOLDER), recursive = TRUE), FALSE)
  
  #MCMC (#RENAME ssir_params AS PARAMS)
  for(i in 2:n_mcmc) {
    
    #print(paste0('i = ', i))
    if(i%%5000 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    ssir_params_dash = c(ssir_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssir_params_dash - vec_min) >= 0){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSIR(epidemic_data, infect_curve_ga, ssir_params_dash, eta)
      #ACCEPTANCE RATIO
      log_accept_ratio = logl_new - log_like
      
      #PRIORS
      log_accept_ratio = log_accept_ratio +
        dexp(ssir_params_dash[1], rate = priors_list$R0_prior[1], log = TRUE) -
        dexp(ssir_params[1], rate = priors_list$R0_prior[1], log = TRUE) +
        dexp(ssir_params_dash[2], rate = priors_list$k_prior[1], log = TRUE) -
        dexp(ssir_params[2], rate = priors_list$k_prior[1], log = TRUE) 
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        ssir_params <- ssir_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssir_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssir_params)
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
      logl_new = LOG_LIKE_SSIR(epidemic_data, infect_curve_ga, ssir_params, eta_dash)
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
      
      ssir_params_matrix[i_thin,] = ssir_params
      log_like_vec[i_thin] <- log_like
      scaling_vec[i_thin] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      eta_matrix[i_thin, ] <- eta 
      sigma_eta_matrix[i_thin, ] = sigma_eta
      i_thin = i_thin + 1
    }
    
  } #END FOR LOOP
  
  #SAVE 
  #saveRDS(epidemic_data, file = paste0(OUTER_FOLDER, 'epidemic_data', seed_count, '.rds' ))
  #saveRDS(ssir_params_matrix, file = paste0(OUTER_FOLDER, 'ssir_params_matrix_', seed_count, '.rds' ))
  #saveRDS(eta_matrix, file = paste0(OUTER_FOLDER, 'eta_matrix_', seed_count, '.rds' ))
  #saveRDS(log_like_vec, file = paste0(OUTER_FOLDER, 'log_like_vec_', seed_count, '.rds' ))
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*num_days)
  
  #Return a, acceptance rate
  return(list(ssir_params_matrix = ssir_params_matrix, eta_matrix = eta_matrix,
              sigma_eta_matrix = sigma_eta_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate, accept_rate_da = accept_rate_da))
} 