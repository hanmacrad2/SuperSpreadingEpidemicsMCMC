#SSEC model

#SIMULATE
SIMULATE_EPI_SSEC <- function(num_days = 50, R0 = 1.2, k = 0.16,
                              shape_gamma = 6, scale_gamma = 1) {
  
  'Simulate an epidemic with Superspreading events
  alpha - R0'

  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    x[t] = rnbinom(1, size = k, mu =  R0*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  x
}

#************************
#* LOG LIKELIHOOD SSEC
#* ***********************
LOG_LIKE_SSEC <- function(x, lambda_vec, ssec_params){
  
  #Params
  R0 = ssec_params[1]; k = ssec_params[2]
  num_days = length(x); loglike = 0
  
  for (t in 2:num_days) {

    loglike = loglike + dnbinom(1, size = k, mu =  R0*lambda_vec[t], log = TRUE)
    
  }
  
  return(loglike)
}

#********************************************************
#1. MCMC INFERENCE FOR SSIC MODEL - INDIVIDUAL R0  (INC. ADAPTIVE SCALING)                           
#********************************************************
MCMC_INFER_SSEC <- function(epidemic_data, n_mcmc,
                            mcmc_inputs = list(mod_start_points = list(m1 = 0.8, m2 = 0.1),
                                               dim = 2, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                               thinning_factor = 10),
                            priors_list = list(r0_prior = c(1,0), k_prior = c(1, 0)),
                            FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper); #NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  num_days = length(epidemic_data)
  vec_min = rep(0, mcmc_inputs$dim)
  count_accept = 0
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #MODEL PARAMETERS
  lambda_vec = get_lambda(epidemic_data)
  ssec_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   #Changed from 0 to NA (As should be overwriting all cases)
  ssec_params_matrix[1,] <- mcmc_inputs$mod_start_points; ssec_params = ssec_params_matrix[1,] #2x1 #as.matrix
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSEC(epidemic_data, scaling_vec, ssec_params);  log_like = log_like_vec[1]
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  scaling_vec <- vector('numeric', mcmc_vec_size); scaling_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssec_params + ssec_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssec_params_matrix[1,]) + tcrossprod(ssec_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #MCMC
  for(i in 2:n_mcmc) {

    if(i%%100 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    ssec_params_dash = c(ssec_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssec_params_dash - vec_min) >= 0){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSEC(epidemic_data, lambda_vec, ssec_params_dash)
      #ACCEPTANCE RATIO
      log_accept_ratio = logl_new - log_like
      
      #PRIORS
      log_accept_ratio = log_accept_ratio -
        priors_list$r0_prior[1]*ssec_params_dash[1] + priors_list$r0_prior[1]*ssec_params[1] - 
        priors_list$k_prior[1]*ssec_params_dash[2] + priors_list$k_prior[1]*ssec_params[2] 
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        ssec_params <- ssec_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssec_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssec_params)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
      
      #ACCEPTANCE PROB - ADAPTIVE SCALING
      accept_prob = min(1, exp(log_accept_ratio))
      
    } else {
     accept_prob = 0
    }
    
    #LAMBDA ADAPTIVE SCALING
    scaling =  scaling*exp(delta/i*(accept_prob - mcmc_inputs$target_acceptance_rate))
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      ssec_params_matrix[i/thinning_factor,] = ssec_params
      log_like_vec[i/thinning_factor] <- log_like
      scaling_vec[i/thinning_factor] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
    }
    
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  
  #SAVE 
  #saveRDS(ssec_params_matrix, file = paste0(OUTER_FOLDER, 'ssec_params_matrix_', seed_count, '.rds' ))
  
  #Return a, acceptance rate
  return(list(ssec_params_matrix = ssec_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate))
} 