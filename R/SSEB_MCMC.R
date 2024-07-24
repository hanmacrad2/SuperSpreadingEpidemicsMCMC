#*************************************
#SSEB MODEL - POISSON-POISSON COMPOUND
#*************************************
#**************************************
#SIMULATE AN EPIDEMIC FROM THE SSEB MODEL
SIMULATE_EPI_SSEB <- function(num_days = 50, r0 = 2.0, alpha = 0.5, beta = 10,
                                 shape_gamma = 6, scale_gamma = 1,
                                 epi_data = c(0,0,0), SIM_DATA = TRUE) {
  
  'Simulate an epidemic with Superspreading events
  alpha - rate of non super-spreading events/days
  Beta = Proportion of superspreading events/days
  beta = increased rate of superspreading event'
  
  #MODEL PARAMS
  gamma = r0*(1 - alpha)/beta #rate of infections = r0_sse/num infections
  #alpha = alpha*r0 #alpha is now the rate
  
  print(paste0('alpha = ', alpha)); print(paste0('beta = ', beta))
  print(paste0('beta = ', beta))
  #Set up
  total_infections = vector('numeric', num_days)
  nsse_infections = vector('numeric', num_days)
  sse_infections = vector('numeric', num_days)
  
  if (SIM_DATA){
    total_infections[1] = 2
    nsse_infections[1] = 2
    sse_infections[1] = 0 
  } else {
    total_infections[1] = epi_data[1]
    nsse_infections[1] = epi_data[1]
    sse_infections[1] = 0 
  }
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infections (tot_rate = lambda) fix notation
    lambda_t = sum(total_infections[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    tot_rate = alpha*r0*lambda_t #Product of infections & their probablilty of infection along the gamma dist at that point in time
    nsse_infections[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, gamma*lambda_t) #Number of super-spreading events (gamma)
    sse_infections[t] = rpois(1, beta*n_t) #z_t: Total infections due to super-spreading event - num of events x Num individuals
    
    total_infections[t] = nsse_infections[t] + sse_infections[t]
  }
  
  total_infections
}

#1. LOG LIKELIHOOD
LOG_LIKE_SSEB <- function(x, lambda_vec, alpha, r0, beta){
  
  #Params
  num_days = length(x); logl = 0
  gamma = r0*(1 - alpha)/beta #r0 = alpha*r0 + beta*gamma (alpha is the proportion of non ss)
  
  for (t in 2:num_days) {
    
    logl = logl + PROBABILITY_XT(x[t], lambda_vec[t], gamma, beta, alpha, r0)
    
  }
  
  return(logl)
}


PROBABILITY_XT <- function(xt, lambda_t, gamma, beta, alpha, r0, max_et = 100) {
  
  'Compound Poisson Prob for SSEB model'
  
  #MAX EVENTS (ADD)
  max_et = GET_MAX_SS_EVENTS(xt) #ADDED 12/01/24
  prob_xt = 0
  vec_prob_xt = vector('numeric', length = max_et + 1)
  
  gamma_lambda_t = gamma*lambda_t
  alpha_r0_lambda_t = alpha*r0*lambda_t
  
  et_values <- 0:max_et
  vec_prob_xt = dpois(et_values, gamma_lambda_t, log = TRUE) +
    dpois(xt, beta*et_values + alpha_r0_lambda_t, log = TRUE)
  
  prob_xt = LOG_SUM_EXP(vec_prob_xt)
  
  return(prob_xt)
}

PROBABILITY_XT_V0 <- function(xt, lambda_t, beta, gamma, alpha, r0, max_et = 100){
  
  'Probability of xt in SSEB model'
  
  #SETUP
  prob_xt = 0
  vec_prob_xt = vector('numeric', length = max_et + 1)
  
  for (et in 0:max_et){
    vec_prob_xt[et + 1] =  dpois(et, beta*lambda_t, log = TRUE) +
      dpois(xt, gamma*et+alpha*r0*lambda_t, log = TRUE)
  }
  
  prob_xt = LOG_SUM_EXP(vec_prob_xt)
  
  return(prob_xt)
}

#GET MAXIMUM POSSIBLE NUMBER OF SSE EVENTS
GET_MAX_SS_EVENTS <- function(data_t) {
  
  if(data_t <= 1000){
    max_et = 100
  } else if (data_t > 1000 & data_t < 10000){
    max_et = 1000
  } else if (data_t > 10000){
    max_et = 10000
  }
}

#PRIOR
SET_SSEB_PRIOR <- function(param, param_dash,
                           alpha_flag = FALSE, r0_flag = FALSE,
                           beta_flag = FALSE){
  
  #PARAMS
  list_priors = GET_LIST_PRIORS_SSEB() 
  PRIORS_USED =  GET_PRIORS_USED() 
  
  if(alpha_flag){
    
    #BETA PRIOR ON a (Just a)
    if (PRIORS_USED$SSEB$alpha$BETA) {
      shape1 = list_priors$alpha[1]
      shape2 = list_priors$alpha[2]
      prior = dbeta(param_dash, shape1, shape2, log = TRUE) -
        dbeta(param, shape1, shape2, log = TRUE) 
    }
    
  } else if (r0_flag) {
    
    if (PRIORS_USED$SSEB$r0$EXP) {
      prior = dexp(param_dash, rate = list_priors$r0[1], log = TRUE) -
        dexp(param, rate = list_priors$r0[1], log = TRUE) 
      
      
    } else if (PRIORS_USED$SSEB$r0$UNIF){
      
      prior = dunif(param_dash, min = list_priors$r0_unif[1], max = list_priors$r0_unif[2], log = TRUE) -
        dunif(param, min = list_priors$r0_unif[1], max = list_priors$r0_unif[2], log = TRUE) 
    }
    
    
  } else if (beta_flag) {
    
    #gamma PRIOR ON beta
    if (PRIORS_USED$SSEB$beta$GAMMA) {
      shape = list_priors$beta[1]
      scale = list_priors$beta[2]
      prior = dgamma(param_dash-1, shape = shape, scale = scale, log = TRUE) -
        dgamma(param-1, shape = shape, scale = scale, log = TRUE) 
    }
  }
  
  return(prior)  
}

#************************************************************************
#1. SSEB MCMC
#************************************************************************
MCMC_INFER_SSEB <- function(epidemic_data, n_mcmc,
                            PRIORS_USED = GET_PRIORS_USED(),
                            param_starts = list(alpha_start = 0.5, beta_start = 10, r0_start = 1.0),
                            mcmc_inputs =  list(alpha_star = 0.3, thinning_factor = 10, burn_in_pc = 0.2), 
                            sigma_starts = list(sigma_alpha = 0.3, sigma_r0 = 0.03, sigma_beta = 3),
                            FLAGS_LIST = list(ADAPTIVE = TRUE, THIN = TRUE, BURN_IN = TRUE)) {
  
  'Returns MCMC samples of SSEB model parameters (alpha, beta, gamma, r0 = alpha + beta*gamma)
  w/ acceptance rate
  Priors
  p(gamma) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #MCMC INITIAL POINTS
  r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  param_starts$r0_start = r0_start
  
  #beta_start = GET_BETA_INITIAL_MCMC(epidemic_data) ?? max num of events etc
  #param_starts$beta_start = beta_start
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  i_thin = 1
  alpha_star = mcmc_inputs$alpha_star
  
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
  alpha_vec[1] <- param_starts$alpha_start
  beta_vec[1] <- param_starts$beta_start
  r0_vec[1] <- r0_start
  gamma_vec[1] <-  r0_vec[1]*(1 - alpha_vec[1])/beta_vec[1]
  log_like_vec[1] <- LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_vec[1], r0_vec[1], beta_vec[1])
  
  #INITIALISE RUNNING PARAMS
  alpha = alpha_vec[1]; beta = beta_vec[1]; 
  gamma = gamma_vec[1]; log_like = log_like_vec[1]
  r0  = r0_vec[1]
  
  #SIGMA
  sigma_alpha = sigma_starts$sigma_alpha; sigma_r0 = sigma_starts$sigma_r0; 
  sigma_beta = sigma_starts$sigma_beta;
  
  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma_alpha_vec <- vector('numeric', mcmc_vec_size); sigma_r0_vec <- vector('numeric', mcmc_vec_size)
    sigma_beta_vec <- vector('numeric', mcmc_vec_size); 
    sigma_alpha_vec[1] =  sigma_alpha; sigma_r0_vec[1] = sigma_r0; sigma_beta_vec[1] =  sigma_beta
    sigma_list = list(sigma_alpha_vec = sigma_alpha_vec, sigma_r0_vec = sigma_r0_vec, 
                      sigma_beta_vec = sigma_beta_vec)
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
    
  } else {
    sigma_list = list(sigma_alpha = sigma_alpha, sigma_r0 = sigma_r0,
                      sigma_beta = sigma_beta)
  }
  
  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0,
                            count_accept_r0 = 0, count_accept3= 0)
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    if (i%%500 == 0) {
      
      print(paste0('i = ', i))
    }
    
    #****************************************************** 
    #alpha
    alpha_dash <- alpha + rnorm(1, sd = sigma_alpha)
    alpha_dash <- alpha_dash %% 2
    if (alpha_dash > 1) {
      alpha_dash <- 2 - alpha_dash
    }
    
    #log a
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_dash, r0, beta)
    
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    
    #Priors
    log_accept_ratio = log_accept_ratio + SET_SSEB_PRIOR(alpha, alpha_dash, alpha_flag = TRUE)
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      alpha <- alpha_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    }
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma_alpha =  sigma_alpha*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************************************************ Only if (b > 0)
    #r0
    r0_dash <- abs(r0 + rnorm(1, sd = sigma_r0))
    
    #loglikelihood
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha, r0_dash, beta)
    log_accept_ratio = logl_new - log_like
    
    #Prior
    log_accept_ratio = log_accept_ratio + SET_SSEB_PRIOR(r0, r0_dash, r0_flag = TRUE)
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      r0 <- r0_dash
      log_like = logl_new
      list_accept_counts$count_accept_r0 = list_accept_counts$count_accept_r0 + 1
      
    }
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma_r0 =  sigma_r0*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #************************************************************************
    #beta
    beta_dash <- beta + rnorm(1, sd = sigma_beta)
    
    if(beta_dash < 1){
      beta_dash = 2 - beta_dash #Prior on c: > 1
    }
    
    #Acceptance Probability
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha, r0, beta_dash)
    log_accept_ratio = logl_new - log_like
    
    #Prior
    log_accept_ratio = log_accept_ratio + SET_SSEB_PRIOR(beta, beta_dash, beta_flag = TRUE)
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      beta <- beta_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma_beta =  sigma_beta*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }
    
    #*****************************************************
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0 && i >= burn_in_start && i_thin <= mcmc_vec_size) {
      alpha_vec[i_thin] <- alpha; r0_vec[i_thin] <- r0
      beta_vec[i_thin] <- beta; gamma_vec[i_thin] <- r0*(1-alpha)/beta
      log_like_vec[i_thin] <- log_like
      sigma_list$sigma_alpha_vec[i_thin] = sigma_alpha; sigma_list$sigma_r0_vec[i_thin] = sigma_r0
      sigma_list$sigma_beta_vec[i_thin] = sigma_beta
      i_thin = i_thin + 1
    }
  }
  
  #Final stats
  accept_rate_alpha = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate_r0 = 100*list_accept_counts$count_accept_r0/(n_mcmc-1) #(list_accept_counts$count_accept_beta + list_reject_counts$count_accept_beta)
  accept_rate_beta = 100*list_accept_counts$count_accept3/(n_mcmc-1)
  
  #Acceptance rates
  list_accept_rates = list(accept_rate_alpha = accept_rate_alpha,
                           accept_rate_r0 = accept_rate_r0, accept_rate_beta = accept_rate_beta)
  print(list_accept_rates)
  
  #OUTPUT
  mcmc_output = list(alpha_vec = alpha_vec, beta_vec = beta_vec, gamma_vec = gamma_vec, r0_vec = r0_vec,
                     log_like_vec = log_like_vec, sigma_list = sigma_list,
                     list_accept_rates = list_accept_rates, 
                     r0_start = r0_start)
  #saveRDS(mcmc_output, file = 'mcmc_sse_output_poisson_compound.rds')
  
  return(mcmc_output)
}
