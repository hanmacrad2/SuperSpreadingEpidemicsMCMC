#********************************************************
#*
#1. SSIB MCMC JOINT
# 
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

#************************
#DATA SIMULATION FUNCTIONS
SIMULATE_EPI_SSIB_LIST = function(num_days = 50, r0 = 2.0,
                                  a = 0.5, b = 10, data_start = c(3,2,1),
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

#***************
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

#****************
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

SET_SSIB_DATA <- function(epidemic_data, FLAGS_DATA_TYPE, data_start = c(3,2,1)){

  'Initialise Super-spreaders & non super-spreaders data'
  
  if(FLAGS_DATA_TYPE$SIM_DATA){
    
    #DATA - AUGMENTED 
    ss = ifelse(epidemic_data > 1, pmax(1, round(0.1*epidemic_data)), 0) #10% of ssib data if > 1, else 1, else 0
    non_ss = epidemic_data - ss
    
    #Starting first time point from 
    non_ss[1] = data_start[2] #Pass in 3,2,1 to first data point
    ss[1] = data_start[3]
    
  } else if (FLAGS_DATA_TYPE$REAL_DATA){
    
    #DATA - AUGMENTED 
    ss = ifelse(epidemic_data > 1, pmax(1, round(0.1*epidemic_data)), 0) #10% of ssib data if > 1, else 1, else 0
    non_ss = epidemic_data - ss
    # 
    # if(epidemic_data[1] == 2){
    #   ss[1] = 0
    #   non_ss[1] = 2
    # }
    
  }
  
  #SS & NS DATA
  print('ss: '); print(ss); print('non_ss: '); print(non_ss)
  output = list(ss = ss, non_ss = non_ss)
  
  return(output)
}

#*************************
#* MCMC 

#' @export
MCMC_INFER_SSIB <- function(epidemic_data, n_mcmc,
                            FLAGS_DATA_TYPE = list(SIM_DATA = TRUE, REAL_DATA = FALSE),
                            PRIORS_USED = GET_PRIORS_USED(),
                           data_start = c(3,2,1),
                            rep_da = 500, mcmc_inputs = list(dim = 3, target_acceptance_rate = 0.234, #0.4 
                                               v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                               thinning_factor = 10, burn_in_pc = 0.2)){    
  
  #INITIALIASE
  param_starts = c(0,0.5,0)
  r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  param_starts[1] = r0_start
  b_start = round(runif(1, 5, 15), 0)
  param_starts[3] = b_start
  
  time = length(epidemic_data) 
  vec_min = c(0, 0, 1); count_accept = 0; 
  vec_accept = vector('numeric', n_mcmc); vec_accept[1] = 0
  
  #THINNING FACTOR + BURN-IN
  thinning_factor = mcmc_inputs$thinning_factor;  i_thin = 1
  mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; print(paste0('N burn-in = ', burn_in_start))
  mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size; print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
  
  #MODEL PARAMETERS
  ssib_params_matrix = matrix(NA, mcmc_vec_size, mcmc_inputs$dim);   
  ssib_params_matrix[1,] <- param_starts; ssib_params = ssib_params_matrix[1,] 
  print(paste0('ssib_params', ssib_params))
  
  #DATA
  print(FLAGS_DATA_TYPE)
  data_ss = SET_SSIB_DATA(epidemic_data, FLAGS_DATA_TYPE, data_start = data_start)
  non_ss = unlist(data_ss$non_ss); ss = unlist(data_ss$ss)
  data = list(non_ss, ss)
  #OUTPUT FROM FUNCTION
  ss_start = ss; ns_start = non_ss 
  
  #STORAGE 
  ns_mcmc = matrix(0, mcmc_vec_size, time)
  ss_mcmc = matrix(0, mcmc_vec_size, time) 
  
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOG_LIKE_SSIB_JOINT(data, ssib_params) #, FLAG_NEGBIN_PARAMATERISATION)
  log_like = log_like_vec[1]
  print(paste0('log_like no1: ', log_like))

  #ADAPTIVE SHAPING PARAMS + VECTORS
  c_star = (2.38^2)/mcmc_inputs$dim; 
  termX = mcmc_inputs$v0 + mcmc_inputs$dim
  delta = 1/(mcmc_inputs$target_acceptance_rate*(1 - mcmc_inputs$target_acceptance_rate))
  x_bar = 0.5*(ssib_params + ssib_params_matrix[1,])
  sigma_i = diag(mcmc_inputs$dim); scaling = 1
  sigma_i = (1/(termX + 3))*(tcrossprod(ssib_params_matrix[1,]) + tcrossprod(ssib_params) -
                               2*tcrossprod(x_bar) + (termX + 1)*sigma_i) #CHANGE TO USE FUNCTIONS
  #Store
  scaling_vec <- vector('numeric', mcmc_vec_size)
  scaling_vec[1] <- 1
  
  #MCMC
  for(i in 2:n_mcmc) {
    
    if(i%%1000 == 0) print(paste0('i = ', i))
    
    #PARAMETERS PROPOSAL
    ssib_params_dash = c(ssib_params +
                           mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    
    #POSTIVE ONLY
    if (min(ssib_params_dash - vec_min) >= 0 && ssib_params_dash[2] < 1){ 
      
      
      #************************************
      #DATA AUGMENTATION
      #***********************************
      #FOR EACH S_T
      for (myrep in 1:rep_da){
        
        #PARAMETERS
        r0 = ssib_params[1];  a = ssib_params[2]; b = ssib_params[3]
        r0_dash = ssib_params_dash[1];  a_dash = ssib_params_dash[2]; b_dash = ssib_params_dash[3]
        
        #PROBABILITIES
        prob = (1 - a)/(a*b + 1 - a)
        prob2 = (1 - a_dash)/(a_dash*b_dash + 1 - a_dash)
        
        #DATA
        I_t = data[[1]] + data[[2]] #'Total epidemic data
        data_dash = data
        
        #PROPOSAL: SS
        q = (data[[2]]/pmax(I_t,1) + prob2)/2  #FUNCTION OF PROB_2; NEW PARAMS. #q = lambda*(s/x + p)/(lambda +1) for lambda = 1 
        s_dash = rbinom(length(I_t), size = I_t, prob = q) #Binomial distribution s_dash
        q2 = (s_dash/pmax(I_t,1) + prob)/2 ##q = lambda*(s/x + p)/(lambda +1) for lambda = 1
        
        data_dash[[2]] = s_dash
        data_dash[[1]] =  I_t - s_dash #n_t = I_t - s_t
        
        #Datapoint: t=1; set to the true/original value
        data_dash[[1]][1]= data[[1]][1]
        data_dash[[2]][1]= data[[2]][1]
        
        #********
        #LOG-LIKELIHOOD -> ACCEPTANCE RATIO
        logl_new = LOG_LIKE_SSIB_JOINT(data_dash, ssib_params_dash)
        
        #ACCEPTANCE RATIO
        log_accept_ratio = logl_new - log_like
        
        #PARAMETER PRIORS
        log_accept_ratio = log_accept_ratio + SET_SSIB_PRIOR_JOINT(ssib_params, ssib_params_dash, PRIORS_USED)
        
        #PROPOSAL RATIOS
        log_accept_ratio = log_accept_ratio + sum(dbinom(data[[2]], size = I_t, prob = q2, log = TRUE)) -
          sum(dbinom(s_dash, size = I_t, prob = q, log = TRUE))   
        
        #METROPOLIS ACCEPTANCE STEP
        if(log(runif(1)) < log_accept_ratio) {
          ssib_params <- ssib_params_dash #PARAMETER UPDATE 
          data <- data_dash #DATA UPDATE 
          log_like <- logl_new
          vec_accept[i] = vec_accept[i] + 1  
        }
        
      }#END DATA AUGMENTAION REPS
      
      #SIGMA - ADAPTIVE SHAPING
      xbar_prev = x_bar
      x_bar = (i-1)/i*xbar_prev + (1/i)*ssib_params
      sigma_i = (1/(i + termX + 1))*( (i + termX)*sigma_i +tcrossprod(ssib_params)
                                      + (i-1)*tcrossprod(xbar_prev)
                                      -i*tcrossprod(x_bar))
      
      #ACCEPTANCE PROB - ADAPTIVE SCALING
      accept_prob = min(1, vec_accept[i]) #1 if accepted, 0 if not accepeted within 50 iterations 
      #accept_prob = min(1, exp(log_accept_ratio))
      
    } else {
      accept_prob = 0
    } #END IF PARAMETER PROPOSAL POSTIVE 
    
    #ADAPTIVE SCALING (needs acceptance probability) #if accept_prob = 1, this increases
    scaling =  scaling*exp(delta/i*(accept_prob - mcmc_inputs$target_acceptance_rate))
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
      ssib_params_matrix[i_thin,] = ssib_params
      log_like_vec[i_thin] <- log_like
      scaling_vec[i_thin] <- scaling #Taking role of sigma, overall scaling constant. Sigma becomes estimate of the covariance matrix of the posterior
      ns_mcmc[i_thin, ] = data[[1]] #Updating all data at once
      ss_mcmc[i_thin, ] = data[[2]]
      i_thin = i_thin + 1
    }
    
  } #END MCMC 
  
  #Final stats
  accept_count = pmin(1, vec_accept) #1 if accepted in that time step at all, 0 otherwise 
  
  
  accept_rate = 100*sum(accept_count)/(n_mcmc-1)
  #accept_da = (100*accept_da)/((n_mcmc-1)*rep_da)
  
  #SS DATA
  ns_median = colMedians(ns_mcmc)
  ss_median = colMedians(ss_mcmc)
  ns_mean = round(colMeans(ns_mcmc))
  ss_mean = round(colMeans(ss_mcmc))
  
  #Return a, acceptance rate
  return(list(ssib_params_matrix = ssib_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate, 
              data_final = data, 
              ns_final = data[[1]], ss_final = data[[2]], #FINAL
              ns_median = ns_median, ss_median = ss_median,
              ns_mean = ns_mean, ss_mean = ss_mean, 
              ss_start = ss_start, ns_start = ns_start,
              ns_mcmc = ns_mcmc, ss_mcmc = ss_mcmc,
              r0_start = r0_start))
} 
