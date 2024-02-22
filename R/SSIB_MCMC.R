#********************************************************
#*
#1. SSIB MCMC JOINT
#* 
#********************************************************

SIMULATE_EPI_SSIB = function(num_days = 50, r0 = 2.0, a = 0.5, b = 10,
                                data_start = c(3,2,1),
                             shape_gamma = 6, scale_gamma = 1) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Params
  c = (r0*(1 - a))/b #r0 = a_prop*r0 + b*c
  #alpha = alpha_prop*r0 
  
  #Set up
  total_infections = vector('numeric', num_days)
  non_ss = vector('numeric', num_days)
  ss = vector('numeric', num_days)
  total_infections[1] = data_start[1] #3
  non_ss[1] =  data_start[2] #2
  ss[1] =  data_start[3] #1 
  
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((non_ss[1:(t-1)] + b*ss[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    non_ss[t] = rpois(1, a*r0*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss[t] = rpois(1, c*lambda_t)
    total_infections[t] = non_ss[t] + ss[t]
  }
  
  total_infections
}

#DATA SIMULATION FUNCTIONS
SIMULATE_EPI_SSIB_LIST = function(num_days = 50, r0 = 2.0,
                                  a = 0.5, b = 10, data_start = c(3,3,1),
                                  shape_gamma = 6, scale_gamma = 1) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Params
  c = (r0*(1 - a))/b #r0 = a_prop*r0 + b*c
  
  #Set up
  total_infections = vector('numeric', num_days)
  non_ss = vector('numeric', num_days)
  ss = vector('numeric', num_days)
  total_infections[1] = data_start[1] #3
  non_ss[1] =  data_start[2] #2
  ss[1] =  data_start[3] #1 
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - 
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((non_ss[1:(t-1)] + b*ss[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    non_ss[t] = rpois(1, a*r0*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss[t] = rpois(1, c*lambda_t)
    total_infections[t] = non_ss[t] + ss[t]
  }
  
  list_ssib_data = list(total_infections = total_infections, non_ss = non_ss,
                        ss = ss)
  
  return(list_ssib_data)
}

#LOG LIKE
LOG_LIKE_SSIB_JOINT <- function(data, ssib_params, 
                                shape_gamma = 6, scale_gamma = 1){
  
  #PARAMS
  r0 = ssib_params[1];  a = ssib_params[2]; b = ssib_params[3]
  
  non_ss = data[[1]]; ss = data[[2]]
  c = (r0*(1 - a))/b
  
  #Params
  num_days = length(non_ss)
  loglike = 0
  
  #INFECTIOUSNESS  - Difference of two GAMMA distributions. Discretized
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) { 
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS
    lambda_t = sum((non_ss[1:(t-1)] + b*ss[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    loglike = loglike + dpois(non_ss[t], a*r0*lambda_t, log = TRUE) +
      dpois(ss[t], c*lambda_t, log = TRUE) 
  }
  
  # if(is.nan(loglike)){
  #   browser()
  # }
  
  return(loglike)
}

#PRIOR
SET_SSIB_PRIOR_JOINT <- function(ssib_params, ssib_params_dash, PRIORS_USED){
  
  #PARAMS
  r0 = ssib_params[1]; r0_dash = ssib_params_dash[1]
  a =  ssib_params[2]; a_dash = ssib_params_dash[2]
  b =  ssib_params[3]; b_dash = ssib_params_dash[3]
  list_priors = GET_LIST_PRIORS_SSIB(); prior = 0
  
  #BETA PRIOR ON a
  if (PRIORS_USED$SSIB$a$BETA) {
    shape1 = list_priors$a[1]
    shape2 = list_priors$a[2]
    prior = dbeta(a_dash, shape1, shape2, log = TRUE) -
      dbeta(a, shape1, shape2, log = TRUE) 
  }
  
  if (PRIORS_USED$SSIB$r0$EXP) {
    prior = prior + dexp(r0_dash, log = TRUE) - dexp(r0, log = TRUE) 
  }
  
  #GAMMA PRIOR ON c
  if (PRIORS_USED$SSIB$b$GAMMA) {
    shape = list_priors$b[1]
    scale = list_priors$b[2]
    prior = prior + dgamma(b_dash -1, shape = shape, scale = scale, log = TRUE) -
      dgamma(b -1, shape = shape, scale = scale, log = TRUE) 
  }
  
  return(prior)  
}

#' @export
MCMC_INFER_SSIB <- function(epidemic_data, n_mcmc, PRIORS_USED = GET_PRIORS_USED(), #list_ssib_data
                                  param_starts = c(2.0, 0.5, 10), data_start = c(3,2,1),
                                  mcmc_inputs = list(dim = 3, target_acceptance_rate = 0.234, #0.4 
                                                     v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                                     thinning_factor = 10, burn_in_pc = 0.2)){    
  
  #INITIALIASE
  r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  param_starts[1] = r0_start
  beta_start = round(runif(1, 5, 15), 0)
  param_starts[3] = beta_start
  
  time = length(epidemic_data) 
  vec_min = c(0, 0, 1); count_accept = 0; 
  
  #THINNING FACTOR + BURN-IN
  thinning_factor = mcmc_inputs$thinning_factor;  i_thin = 1
  mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; print(paste0('N burn-in = ', burn_in_start))
  mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size; print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
  
  #MODEL PARAMETERS
  ssib_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   
  ssib_params_matrix[1,] <- param_starts; ssib_params = ssib_params_matrix[1,] 
  print(paste0('ssib_params', ssib_params))
  
  #DATA - AUGMENTED 
  ss = ifelse(epidemic_data > 1, pmax(1, round(0.1*epidemic_data)), 0)
  non_ss = epidemic_data - ss
  #Starting first time point from the truth
  non_ss[1] = data_start[2]
  ss[1] = data_start[3]
  ss_start = ss; ns_start = non_ss
  data = list(non_ss, ss)
  #data = list(unlist(list_ssib_data$non_ss), unlist(list_ssib_data$ss))
  print('ss: '); print(ss); print('non_ss: '); print(non_ss)
  
  #SSI - DATA AUGMENTATION STORE
  non_ss = matrix(0, mcmc_vec_size, time)
  ss = matrix(0, mcmc_vec_size, time) 
  vec_accept_da = vector('numeric', length = time)
  accept_da = 0
  
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSIB_JOINT(data, ssib_params) #, FLAG_NEGBIN_PARAMATERISATION)
  log_like = log_like_vec[1]
  print(paste0('log_like 0', log_like))
  
  #ADAPTIVE SHAPING PARAMS + VECTORS
  scaling_vec <- vector('numeric', mcmc_vec_size); scaling_vec[1] <- 1
  c_star = (2.38^2)/mcmc_inputs$dim; termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssib_params + ssib_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssib_params_matrix[1,]) + tcrossprod(ssib_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  
  #MCMC
  for(i in 2:n_mcmc) {
    
    if(i%%1000 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    ssib_params_dash = c(ssib_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssib_params_dash - vec_min) >= 0 && ssib_params_dash[2] < 1){ 
      
      #LOG LIKELIHOOD
      logl_new = LOG_LIKE_SSIB_JOINT(data, ssib_params_dash) 
      
      #ACCEPTANCE RATIO
      log_accept_ratio = logl_new - log_like
      
      #PRIORS
      log_accept_ratio = log_accept_ratio + SET_SSIB_PRIOR_JOINT(ssib_params, ssib_params_dash, PRIORS_USED)
      
      #METROPOLIS ACCEPTANCE STEP
      if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
        ssib_params <- ssib_params_dash
        count_accept = count_accept + 1
        log_like = logl_new
      }
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssib_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssib_params)
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
      ssib_params_matrix[i_thin,] = ssib_params
      log_like_vec[i_thin] <- log_like
      scaling_vec[i_thin] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      i_thin = i_thin + 1
    }
    
    #************************************
    #DATA AUGMENTATION
    #***********************************
    #FOR EACH S_T
    for (myrep in 1:10){
      
      data_dash = data
      I_t = data[[1]] + data[[2]] #x = 'epidemic data'
      
      #PARAMS
      r0 = ssib_params[1];  a = ssib_params[2]; b = ssib_params[3]
      prob = (1 - a)/(a*b + 1 - a)
      proposal_ss = rbinom(length(I_t), size = I_t, prob = prob) #Binomial distribution; from the total infections
      
      data_dash[[2]] = proposal_ss
      data_dash[[1]] =  I_t - proposal_ss #n_t = I_t - s_t
      
      #Keep first data == truth
      data_dash[[1]][1]= data[[1]][1]
      data_dash[[2]][1]= data[[2]][1]
      
      logl_new = LOG_LIKE_SSIB_JOINT(data_dash, ssib_params) 
      log_accept_ratio = logl_new - log_like - sum(dbinom(proposal_ss, size = I_t, prob = prob, log = TRUE)) + 
       sum(dbinom(data[[2]], size = I_t, prob = prob, log = TRUE))
        
        #METROPOLIS ACCEPTANCE STEP
        if(log(runif(1)) < log_accept_ratio) {
          data <- data_dash
          log_like <- logl_new
          accept_da = accept_da + 1
          #vec_accept_da[t] =  vec_accept_da[t] + 1
        }
    }
      
      #Store
      if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
        non_ss[i_thin, ] = data[[1]] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i_thin, ] = data[[2]]
      }
    }
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  accept_da = (100*accept_da)/((n_mcmc-1)*10)
  
  #Return a, acceptance rate
  return(list(ssib_params_matrix = ssib_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate, accept_da = accept_da,
              data = data, ss_inf = data[[1]], ns_inf = data[[2]],
              ss_start = ss_start, ns_start = ns_start,
              non_ss_tot = non_ss, ss_tot = ss,
  r0_start = r0_start))
} 
