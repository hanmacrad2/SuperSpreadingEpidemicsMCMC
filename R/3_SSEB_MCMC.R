#*************************************
#SSEB MODEL - POISSON-POISSON COMPOUND
#*************************************

#**************************************
#SIMULATE AN EPIDEMIC FROM THE SSEB MODEL
SIMULATE_EPI_SSEB <- function(num_days = 30, alphaX = 0.8, betaX = 0.2, gammaX = 10,
                              shape_gamma = 6, scale_gamma = 1,
                              epi_data = c(0,0,0), SIM_DATA = TRUE) {
  
  'Simulate an epidemic with Superspreading events
  alpha - rate of non super-spreading events/days
  Beta = Proportion of superspreading events/days
  gamma = increased rate of superspreading event'
  
  #MODEL PARAMS
  print(paste0('alpha = ', alphaX)); print(paste0('beta = ', betaX))
  print(paste0('gamma = ', gammaX))
  #Set up
  total_infecteds = vector('numeric', num_days)
  nsse_infecteds = vector('numeric', num_days)
  sse_infecteds = vector('numeric', num_days)
  
  if (SIM_DATA){
    total_infecteds[1] = 2
    nsse_infecteds[1] = 2
    sse_infecteds[1] = 0 
  } else {
    total_infecteds[1] = epi_data[1]
    nsse_infecteds[1] = epi_data[1]
    sse_infecteds[1] = 0 
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(total_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, betaX*lambda_t) #Number of super-spreading events (beta)
    sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}

#1. LOG LIKELIHOOD
LOG_LIKE_SSEB <- function(x, lambda_vec, alpha_prop, r0, gamma){
  
  #Params
  num_days = length(x); logl = 0
  
  beta = r0*(1 - alpha_prop)/gamma #r0 = alpha*r0 + beta*gamma (alpha is the proportion of non ss)
  
  alpha = alpha_prop*r0 #alpha is now the rate
   
  for (t in 2:num_days) {
    
    #lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum_xt = 0
    term1 = exp(-alpha*lambda_vec[t]); term2 = alpha*lambda_vec[t]
    
    for (nt in 0:x[t]){ #Sum for all possible values of nt, xt-nt
      
      #Log likelihood
      st = x[t] - nt
      inner_sum_xt = inner_sum_xt + 
        dpois(nt, alpha*lambda_vec[t])*
        PROBABILITY_ST(st, lambda_vec[t], alpha, beta, gamma)
    } 
    
    logl = logl + log(inner_sum_xt) 
  }
  
  return(logl)
}

#2. PROBABILITY OF ST
PROBABILITY_ST <- function(st, lambda_t, alphaX, betaX, gammaX, max_et = 5){
  
  'Probability of St'
  prob_st = 0
  
  for (et in 0:max_et){
    prob_st = prob_st + dpois(et, betaX*lambda_t)*dpois(st, gammaX*et)
  }
  
  return(prob_st)
}

#PRIOR
SET_SSEB_PRIOR <- function(param, param_dash, r0_temp,
                           list_priors, PRIORS_USED,
                           alpha_flag = FALSE, r0_flag = FALSE, gamma_flag = FALSE){
  
  if(alpha_flag){
    
    #BETA PRIOR ON a (Just a)
    if (PRIORS_USED$SSEB$alpha$BETA) {
      shape1 = list_priors$alpha[1]
      shape2 = list_priors$alpha[2]
      p = dbeta(param_dash, shape1, shape2, log = TRUE) -
        dbeta(param, shape1, shape2, log = TRUE) 
    }
    
  } else if (r0_flag) {
    
    if (PRIORS_USED$SSEB$R0$EXP) {
      p = dexp(param_dash, log = TRUE) - dexp(param, log = TRUE) 
    }
    
  } else if (gamma_flag) {
    
    #GAMMA PRIOR ON c
    if (PRIORS_USED$SSEB$gamma$GAMMA) {
      shape = list_priors$gamma[1]
      scale = list_priors$gamma[2]
      #browser()
      p = dgamma(param_dash -1, shape = shape, scale = scale, log = TRUE) -
        dgamma(param-1, shape = shape, scale = scale, log = TRUE) 
    }
  }
  
  return(p)  
}

#************************************************************************
#1. SSEB MCMC
#************************************************************************
MCMC_INFER_SSEB <- function(epidemic_data, n_mcmc = 30000,
                            mcmc_inputs = 
                              list(param_starts = list(alpha_start = 0.8, beta_start = 0.1, gamma_start = 10),
                                   alpha_star = 0.4, thinning_factor = 10, burn_in_pc = 0.2), 
                            sigma_starts = list(sigma_alpha = 0.3, sigma_r0 = 0.03, sigma_gamma = 3),
                            FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE, BURN_IN = TRUE)) {
  
  'Returns MCMC samples of SSEB model parameters (alpha, beta, gamma, r0 = alpha + beta*gamma)
  w/ acceptance rate
  Priors
  p(gamma) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  
  #PRIORS
  list_priors = GET_LIST_PRIORS_SSEB() 
  PRIORS_USED =  GET_PRIORS_USED() 
  
  #THINNING FACTOR
  i_thin = 1
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #BURN_IN: Initial samples are not completely valid; the Markov Chain has not stabilized to the stationary distribution. The burn in samples allow you to discard these initial samples that are not yet at the stationary.
  if(FLAGS_LIST$BURN_IN){
    burn_in_start = mcmc_inputs$burn_in_pc*n_mcmc; print(paste0('N burn-in = ', burn_in_start))
    
    mcmc_vec_size = mcmc_vec_size - mcmc_inputs$burn_in_pc*mcmc_vec_size;
    print(paste0('Post burn-in mcmc vec size = ', mcmc_vec_size))
  }
  
  #INITIALISE MCMC VECTORS
  alpha_vec <- vector('numeric', mcmc_vec_size); beta_vec <- vector('numeric', mcmc_vec_size)
  gamma_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);
  
  #INITIALISE MCMC[1]
  alpha_vec[1] <- mcmc_inputs$param_starts$alpha_start; beta_vec[1] <- mcmc_inputs$param_starts$beta_start
  gamma_vec[1] <- mcmc_inputs$param_starts$gamma_start; r0_vec[1] <- alpha_vec[1] + beta_vec[1]*gamma_vec[1]
  log_like_vec[1] <- LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_vec[1], r0_vec[1], gamma_vec[1])
  
  #INITIALISE RUNNING PARAMS
  alpha = alpha_vec[1]; beta = beta_vec[1]; 
  gamma = gamma_vec[1]; log_like = log_like_vec[1]
  r0  = r0_vec[1]
  
  #SIGMA
  sigma_alpha = sigma_starts$sigma_alpha; sigma_r0 = sigma_starts$sigma_r0; 
  sigma_gamma = sigma_starts$sigma_gamma;
  
  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma_alpha_vec <- vector('numeric', mcmc_vec_size); sigma_r0_vec <- vector('numeric', mcmc_vec_size)
    sigma_gamma_vec <- vector('numeric', mcmc_vec_size); 
    sigma_alpha_vec[1] =  sigma_alpha; sigma_r0_vec[1] = sigma_r0; sigma_gamma_vec[1] =  sigma_gamma
    sigma_list = list(sigma_alpha_vec = sigma_alpha_vec, sigma_r0_vec = sigma_r0_vec, 
                      sigma_gamma_vec = sigma_gamma_vec)
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
    
  } else {
    
    sigma_list = list(sigma_alpha = sigma_alpha, sigma_r0 = sigma_r0,
                      sigma_gamma = sigma_gamma)
  }
  
  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0,
                            count_accept_r0 = 0, count_accept3= 0)
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    if (i%%1000 == 0) {
      
      print(paste0('i = ', i))
    }
    
    #****************************************************** 
    #alpha
    alpha_dash <- alpha + rnorm(1, sd = sigma_alpha)
    
    while(alpha_dash < 0 || alpha_dash > 1){
      
      if (alpha_dash > 1){ #keep alpha between 0 and 1 
        alpha_dash = 2 - alpha_dash
      }
      alpha_dash = abs(alpha_dash) 
    }
    
    #log a
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_dash, r0, gamma)
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    
    #Priors
    log_accept_ratio = log_accept_ratio + SET_SSEB_PRIOR(alpha, alpha_dash, r0,
                                                         list_priors, PRIORS_USED, alpha_flag = TRUE)
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      alpha <- alpha_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
      
      #Sigma (Adaptive)
      if (FLAGS_LIST$ADAPTIVE){
        accept_prob = min(1, exp(log_accept_ratio))
        sigma_alpha =  sigma_alpha*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
      }
    }
    
    #************************************************************************ Only if (b > 0)
    #r0
    r0_dash <- r0 + rnorm(1, sd = sigma_r0)
    
    #loglikelihood
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha, r0_dash, gamma)
    log_accept_ratio = logl_new - log_like
    
    #Prior
    r0_temp = r0
    log_accept_ratio = log_accept_ratio + SET_SSEB_PRIOR(r0, r0_dash, r0_temp, 
                                                         list_priors, PRIORS_USED, r0_flag = TRUE)
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      r0 <- r0_dash
      log_like = logl_new
      list_accept_counts$count_accept_r0 = list_accept_counts$count_accept_r0 + 1
      
      #Sigma (Adpative)
      if (FLAGS_LIST$ADAPTIVE){
        accept_prob = min(1, exp(log_accept_ratio))
        sigma_r0 =  sigma_r0*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
      }
    }
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma + rnorm(1, sd = sigma_gamma)
    
    if(gamma_dash < 1){
      
      gamma_dash = 2 - gamma_dash #Prior on c: > 1
    }
    
    #Acceptance Probability
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha, r0, gamma_dash)
    log_accept_ratio = logl_new - log_like
    
    #Prior
    log_accept_ratio = log_accept_ratio + SET_SSEB_PRIOR(gamma, gamma_dash, r0,
                                                         list_priors, PRIORS_USED, gamma_flag = TRUE)
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      gamma <- gamma_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
      
      #Sigma (Adpative)
      if (FLAGS_LIST$ADAPTIVE){
        accept_prob = min(1, exp(log_accept_ratio))
        sigma_gamma =  sigma_gamma*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
      }
    }
    
    #*****************************************************
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
      alpha_vec[i_thin] <- alpha; r0_vec[i_thin] <- r0
      gamma_vec[i_thin] <- gamma; beta_vec[i_thin] <- (r0 - alpha)/gamma
      log_like_vec[i_thin] <- log_like
      sigma_list$sigma_alpha_vec[i_thin] = sigma_alpha; sigma_list$sigma_r0_vec[i_thin] = sigma_r0
      sigma_list$sigma_gamma_vec[i_thin] = sigma_gamma#
      i_thin = i_thin + 1
    }
  }
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate_r0 = 100*list_accept_counts$count_accept_r0/(n_mcmc-1) #(list_accept_counts$count_accept_beta + list_reject_counts$count_accept_beta)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1)
  
  #Acceptance rates
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate_r0 = accept_rate_r0, accept_rate3 = accept_rate3)
  print(list_accept_rates)
  
  #OUTPUT
  mcmc_output = list(alpha_vec = alpha_vec, beta_vec = beta_vec, gamma_vec = gamma_vec, r0_vec = r0_vec,
                     log_like_vec = log_like_vec, sigma_list = sigma_list,
                     list_accept_rates = list_accept_rates)
  #saveRDS(mcmc_output, file = 'mcmc_sse_output_poisson_compound.rds')
  
  return(mcmc_output)
}

