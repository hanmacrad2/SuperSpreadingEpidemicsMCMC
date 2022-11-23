#********************************************************
#1. INDIVIDUAL R0 MCMC ADAPTIVE SHAPING                           
#********************************************************
library(MASS)

#*********************************************
#* SIMULATE ETA
#**********************************************
SIMULATE_SSID = function(num_days = 110, alphaX = 1.2, k = 0.16,
                        shape_gamma = 6, scale_gamma = 1) {
  
  'Simulate from the Negative Binomial model'
  
  #INTIALISE VECTORS
  x = vector('numeric', num_days); x[1] = 2
  eta_vec = vector('numeric', num_days); 
  
  #INFECTIOUSNESS (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #DAYS OF THE EPIDEMIC
  for (t in 2:num_days) {
    
    #ETA (t-1)
    eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = alphaX/k) #Draw eta from previous time step
    #INFECTIVITY
    infectivity = rev(prob_infect[1:t]) 
    #POISSON; OFFSPRINT DISTRIBUTION
    total_rate = sum(eta_vec*infectivity) #DOT PRODUCT
    x[t] = rpois(1, total_rate)
    
  }
  return(list(epidemic_data = x, eta_vec = eta_vec))
}

#LOG LIKELIHOOD
#**********************************************
LOG_LIKELIHOOD_SSID <- function(x, infectivity_vec, ssid_params, eta){ #eta - a vector of length x. eta[1] = infectivity of xt[1]
  
  #Params
  num_days = length(x)
  alpha = ssid_params[1]; k = ssid_params[2]
  loglike = 0; count_inf = 0
  
  for (t in 2:num_days) {
    
    #INFECTIVITY
    #infectivity = rev(prob_infect[1:t-1]) #Current infectivity dependent on people already infected #rev(prob_infect[1:(t-1)]) 
    total_rate = sum(eta[1:(t-1)]*infectivity_vec) 
    
    if(x[t-1]==0){ #No infections 
        if(eta[t-1] == 0) eta_prob = 0 #Don't change likelihood
        else eta_prob = -Inf #Make whole likelihood zero 
    } else { #If there IS infections
      eta_prob = dgamma(eta[t-1], shape = x[t-1]*k, scale = alpha/k, log = TRUE)
    }
    
    loglike = loglike + eta_prob 
    
    poi_prob = x[t]*log(total_rate) - total_rate #- lfactorial(x[t]) CANCELS IN ACCEPTANCE RATIO
    
    if (!is.na(poi_prob)) loglike = loglike + poi_prob #Note want to include -Inf so don't filter infinite values
  }
  
  return(loglike)
}

#********************************************************
#1. INDIVIDUAL R0 MCMC ADAPTIVE SHAPING                           
#********************************************************
MCMC_ADAPTIVE_SSID <- function(dataX, OUTER_FOLDER, seed_count,
                              mcmc_inputs = list(n_mcmc = 20000, #00,
                                                 mod_start_points = c(1.2, 0.16),
                                                 dim = 2, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                                 thinning_factor = 10),
                              priors_list = list(k_prior = c(1, 0)),
                              FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE, PRIOR_K1 = TRUE,
                                                PRIOR_K2 = FALSE)) {    
  
  #NOTE:
  #i - 1 = n (Simon's paper); #NOTE NO REFLECTION, NO TRANSFORMS, MORE INTELLIGENT ADAPTATION
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  
  #MCMC PARAMS + VECTORS
  num_days = length(dataX); n_mcmc = mcmc_inputs$n_mcmc;
  vec_min = rep(0, mcmc_inputs$dim);
  count_accept = 0; count_accept_da = 0
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #MODEL PARAMETERS
  infectivity_vec = get_infectivity(dataX)
  eta = dataX; eta_matrix = matrix(NA, mcmc_vec_size, num_days); 
  ssid_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   #Changed from 0 to NA (As should be overwriting all cases)
  ssid_params_matrix[1,] <- mcmc_inputs$mod_start_points; ssid_params = ssid_params_matrix[1,] #2x1 #as.matrix
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKELIHOOD_SSID(dataX, infectivity_vec, ssid_params, eta);  log_like = log_like_vec[1]
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  lambda_vec <- vector('numeric', mcmc_vec_size); lambda_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssid_params + ssid_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); lambda_i = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssid_params_matrix[1,]) + tcrossprod(ssid_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #SIGMA ETA (ONE DIMENSION - ROBINS MUNROE) 
  sigma_eta =  0.5*rep(1, num_days)
  sigma_eta_matrix = matrix(0, mcmc_vec_size, num_days); sigma_eta_matrix[1,] =  sigma_eta;
  
  #DIRECTORY - SAVING
  ifelse(!dir.exists(file.path(OUTER_FOLDER)), dir.create(file.path(OUTER_FOLDER), recursive = TRUE), FALSE)
  
  #MCMC (#RENAME ssid_params AS PARAMS)
  for(i in 2:n_mcmc) {
    
    #print(paste0('i = ', i))
    if(i%%100 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    ssid_params_dash = c(ssid_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = lambda_i*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssid_params_dash - vec_min) >= 0){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_SSID(dataX, infectivity_vec, ssid_params_dash, eta)
      #ACCEPTANCE RATIO
      log_accept_ratio = logl_new - log_like
      
      #PRIOR (*Currently improper prior on alpha. Zero weight on alpha being less than 1, infinite weight on alpha being > 1)
      #Prior on alpha currently 0; log(0); 1
      if (FLAGS_LIST$PRIOR_K1){
        log_accept_ratio = log_accept_ratio - ssid_params_dash[2] + ssid_params[2]  #exp(1) prior on k 
      } else if (FLAGS_LIST$PRIOR_K2){
        log_accept_ratio = log_accept_ratio - priors_list$k_prior[1]*ssid_params_dash[2] + priors_list$k_prior[1]* ssid_params[2] 
      }
      
      #ALPHA PRIOR
      #SAME AMOUNT OF MASS <1>. 
      #Gamma with median 1 Gamma(2,1) Mode:1. or exp(1)
      log_accept_ratio = log_accept_ratio - ssid_params_dash[1] + ssid_params[1]
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        ssid_params <- ssid_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssid_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssid_params)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
      
      #ACCEPTANCE PROB - LAMBDA ADAPTIVE SCALING
      accept_prob = min(1, exp(log_accept_ratio))
      
    } else {
      
      accept_prob = 0
    }
    
    #LAMBDA ADAPTIVE SCALING
    lambda_i =  lambda_i*exp(delta/i*(accept_prob - mcmc_inputs$target_acceptance_rate))
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    for(t in which(dataX[1:(num_days-1)] > 0)){ #Only update eta's where x[t] > 0
      
      v = rep(0, length(eta)); v[t] = 1
      
      #METROPOLIS STEP 
      eta_dash = abs(eta + rnorm(1,0, sigma_eta[t])*v) #normalise the t_th element of eta #or variance = x[t]
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKELIHOOD_SSID(dataX, infectivity_vec, ssid_params, eta_dash)
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
    if (i%%thinning_factor == 0) {
      ssid_params_matrix[i/thinning_factor,] = ssid_params
      log_like_vec[i/thinning_factor] <- log_like
      lambda_vec[i/thinning_factor] <- lambda_i #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      eta_matrix[i/thinning_factor, ] <- eta 
      sigma_eta_matrix[i/thinning_factor, ] = sigma_eta
    }
    
  } #END FOR LOOP
  
  #SAVE 
  #saveRDS(dataX, file = paste0(OUTER_FOLDER, 'dataX', seed_count, '.rds' ))
  #saveRDS(ssid_params_matrix, file = paste0(OUTER_FOLDER, 'ssid_params_matrix_', seed_count, '.rds' ))
  #saveRDS(eta_matrix, file = paste0(OUTER_FOLDER, 'eta_matrix_', seed_count, '.rds' ))
  #saveRDS(log_like_vec, file = paste0(OUTER_FOLDER, 'log_like_vec_', seed_count, '.rds' ))
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_rate_da = 100*count_accept_da/((n_mcmc-1)*num_days)
  
  #Return a, acceptance rate
  return(list(ssid_params_matrix = ssid_params_matrix, eta_matrix = eta_matrix,
              sigma_eta_matrix = sigma_eta_matrix,
              log_like_vec = log_like_vec, lambda_vec = lambda_vec, 
              accept_rate = accept_rate, accept_rate_da = accept_rate_da))
} 