#********************************************************
#*
#1. SSIB MCMC JOINT
#* 
#********************************************************

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
MCMC_INFER_SSIB_JOINT <- function(epidemic_data, list_ssib_data, n_mcmc, PRIORS_USED = GET_PRIORS_USED(),
                           param_starts = c(3.0, 0.85, 10),
                           mcmc_inputs = list(dim = 3, target_acceptance_rate = 0.4, v0 = 100,  #priors_list = list(alpha_prior = c(1, 0), k_prior = c()),
                                              thinning_factor = 10, burn_in_pc = 0.2)){    
  
  #INITIALIASE
  #r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  #param_starts[1] = r0_start
  r0_start = param_starts[1]
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
  
  #SSI - DATA AUGMENTATION PARAMS
  # ss = ifelse(epidemic_data > 1, pmax(1, round(0.1*epidemic_data)), 0)
  # non_ss = epidemic_data - ss
  #ss_start = ss; ns_start = non_ss
  # 
  # #DATA USED!!
  # data = list(non_ss, ss)
  data = list(unlist(list_ssib_data$non_ss), unlist(list_ssib_data$ss))
  #print('ss: '); print(ss); print('non_ss: '); print(non_ss)
  
  #SSI - DATA AUGMENTATION STORE
  non_ss = matrix(0, mcmc_vec_size, time)
  ss = matrix(0, mcmc_vec_size, time) 
  vec_accept_da = vector('numeric', length = time)
  
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
    for(t in 1:time){
      
      data_dash = data
      
      #Propose
      x = data[[1]][t] + data[[2]][t]
      d = max(1, round(runif(1, min = 0, max = x/10)))
      
      #STOCHASTIC PROPOSAL for s
      if (runif(1) < 0.5) {
        st_dash = data[[2]][t] + d #1
      } else {
        st_dash = data[[2]][t] - d #1
      }
      
      #ACCEPTANCE PROBABILITY
      data_dash[[2]][t] = st_dash #s_t = st_dash
      data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t
      
      #CRITERIA FOR S_T & N_T
      if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
        log_accept_ratio = -Inf 
      } else {
        logl_new = LOG_LIKE_SSIB_JOINT(data_dash, ssib_params) #+include alpha prior 
        log_accept_ratio = logl_new - log_like
      }
      
      #METROPOLIS ACCEPTANCE STEP
      if(log(runif(1)) < log_accept_ratio) {
        data <- data_dash
        log_like <- logl_new
        vec_accept_da[t] =  vec_accept_da[t] + 1
      }
      
      #Store
      if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
        non_ss[i_thin, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i_thin, t] = data[[2]][t]
      }
    }
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  
  #Return a, acceptance rate
  return(list(ssib_params_matrix = ssib_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate, vec_accept_da = vec_accept_da,
              data = data, ss_inf = data[[1]], ns_inf = data[[2]],
              #ss_start = ss_start, ns_start = ns_start,
              non_ss_tot = non_ss, ss_tot = ss,
              r0_start = r0_start))
} 
