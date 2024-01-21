#SSIB ORIG PARAM (epidemic_modelling/model_ssi/1_SSI_model_w_data_aug_version_total.R)

#*******************************************************
#I SUPER-SRPEADING INDIVIDUALS (SUPER-SPREADERS) SIMULATION 
simulation_super_spreaders = function(num_days = 50, a = 0.8, b = 0.1, c = 10, 
                                      shape_gamma = 6, scale_gamma = 1) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 3
  nss_infecteds[1] = 2
  ss_infecteds[1] = 1 
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((total_infecteds[1:(t-1)] + c*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, a*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, b*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  return(list(nss_infecteds, ss_infecteds))
}

#****************************************************************
#II. MODEL SSI - LOG LIKELIHOOD
#****************************************************************
LOG_LIKE_SSI_ORIG <- function(sim_data, a, b, c){
  
  #Data
  n = sim_data[[1]]; s = sim_data[[2]]
  
  #Params
  num_days = length(n)
  shape_gamma = 6; scale_gamma = 1
  logl = 0
  
  #INFECTIOUSNESS  - Difference of 2 GAMMA distributions. Discretized 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 1:num_days) { #*1 or 2
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS 
    lambda_t = sum((n[1:(t-1)] + c*s[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    #LOG-LIKELIHOOD 
    logl = logl - lambda_t*(a + b) + n[t]*(log(a) + log(lambda_t)) + s[t]*(log(b) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
  }
  
  logl
}



#************************************************************************
#III. SSI MCMC                              (W/ DATA AUGMENTATION OPTION)
#************************************************************************

MCMC_SSIB_ORIG <- function(data, n_mcmc, 
                           model_params = c(0.8, 0.1, 10), thinning_factor = 10, burn_in_pc = 0.2,
                     priors_list = list(a_prior = c(1, 0), b_prior = c(10, 1/100), b_prior_exp = c(1,0),
                                        c_prior = c(10, 1), c_prior_exp = c(0.1,0)),
                     FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE, PRIOR = TRUE,
                                       B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE)) { 
  
  'Returns MCMC samples of SSI model parameters (a, b, c, r0 = a + b*c) 
  w/ acceptance rates.
  INCLUDES; DATA AUGMENTATION, B-C transform' 
  print('MCMC SUPERSPREADER INDIVIDUALS')
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS (THINK ABOUT STARTING POINTS)
  #**********************************************
  time = length(data[[1]]); i_thin = 1
  mcmc_vec_size = n_mcmc/thinning_factor #Eg 125k/10 = 12.5k
  mcmc_vec_size = mcmc_vec_size - burn_in_pc*mcmc_vec_size #12.5k-2.5k = 10k :)
  mcmc_vec_size = mcmc_vec_size + 1
  burn_in_start = burn_in_pc*n_mcmc
  print(paste0('mcmc_vec_size: ', mcmc_vec_size))
  #MCMC VECTORS
  a_vec <- vector('numeric', mcmc_vec_size); b_vec <- vector('numeric', mcmc_vec_size)
  c_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size)
  
  #INITIALISE: RUNNING PARAMS
  a = model_params[1]; b =  model_params[2]; 
  c = model_params[3]; log_like = LOG_LIKE_SSI_ORIG(data, a, b, c)
  
  #INITIALISE: MCMC[1] of MCMC VECTORS
  a_vec[1] <- a; b_vec[1] <- b
  c_vec[1] <- c; r0_vec[1] <- a + b*c
  log_like_vec[1] <- log_like
  
  #SIGMA (NOT ADAPTIVE)
  sigma_a = 0.25*a; sigma_b = 0.25*b    
  sigma_c = 0.05*c; sigma_bc = 0.05*c
  sigma = list(sigma_a = sigma_a, sigma_b = sigma_b, #Acc rate too big -> Make sigma bigger. 
               sigma_c = sigma_c, sigma_bc = sigma_bc) #Acc rate too small -> make sigma smaller
  
  #INITIALISE: ACCEPTANCE COUNTS 
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0)
  list_reject_counts = list(count_reject2 = 0, count_reject3 = 0,
                            count_reject4 = 0)
  
  mat_count_da = matrix(0, mcmc_vec_size, time) #i x t
  non_ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  ss = matrix(0, mcmc_vec_size, time) #USE THINNING FACTOR
  vec_accept_da = vector('numeric', length = time)
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2: n_mcmc) {
    
    if(i%%1000 == 0){ print(paste0('i:', i, ' loglike: ', log_like))}
    #****************************************************** 
    #a
    a_dash <- a + rnorm(1, sd = sigma$sigma_a) 
    if(a_dash < 0){
      a_dash = abs(a_dash)
    }
    
    #log a
    logl_new = LOG_LIKE_SSI_ORIG(data, a_dash, b, c)
    log_accept_prob = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_prob = log_accept_prob - a_dash + a
    }
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      a <- a_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    } 
    
    #************************************************************************ Only if (b > 0){ ?
    #b  
    b_dash <- b + rnorm(1, sd = sigma$sigma_b) 
    if(b_dash < 0){
      b_dash = abs(b_dash)
    }
    #loglikelihood
    logl_new = LOG_LIKE_SSI_ORIG(data, a, b_dash, c)
    log_accept_prob = logl_new - log_like
    
    #Priors
    if (FLAGS_LIST$B_PRIOR_GAMMA){
      log_accept_prob = log_accept_prob +
        dgamma(b_dash, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE) -
        dgamma(b, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE)
    } else {
      log_accept_prob = log_accept_prob - b_dash + b 
    }
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      b <- b_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    } 
    
    #************************************************************************
    #c
    c_dash <- c + rnorm(1, sd = sigma$sigma_c) 
    if(c_dash < 1){
      c_dash = 2 - c_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSI_ORIG(data, a, b, c_dash)
    log_accept_prob = logl_new - log_like 
    
    #Priors
    if(FLAGS_LIST$C_PRIOR_GAMMA){
      log_accept_prob = log_accept_prob + dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[1], log = TRUE) -
        dgamma(c, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
    } else {
      log_accept_prob = log_accept_prob - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c
    }
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      c <- c_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 =  list_accept_counts$count_accept3 + 1
    }
    
    #*****************************************************
    #B-C TRANSFORM
    if(FLAGS_LIST$BC_TRANSFORM){
      
      c_dash <- c + rnorm(1, sd = sigma$sigma_bc)
      #Prior > 1
      if(c_dash < 1){
        c_dash = 2 - c_dash
      }
      #New b
      b_transform = ((a + b*c) - a)/c_dash #b = (r0 - a)c
      
      if(b_transform >= 0){ #Only accept values of b > 0
        
        logl_new = LOG_LIKE_SSI_ORIG(data, a, b_transform, c_dash)
        log_accept_prob = logl_new - log_like
        
        #PRIORS
        #b prior
        if (FLAGS_LIST$B_PRIOR_GAMMA) {
          tot_b_prior = dgamma(b_transform, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE) -
            dgamma(b, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE)
        } else { 
          tot_b_prior = - b_transform + b #exp(1) piror
        }
        
        #c prior
        if (FLAGS_LIST$C_PRIOR_GAMMA) {
          tot_c_prior = dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE) -
            dgamma(b, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
        } else { 
          tot_c_prior = priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c 
        }
        
        #LOG ACCEPT PROB
        log_accept_prob = log_accept_prob + tot_b_prior + tot_c_prior 
        
        #Metropolis Step
        if (log(runif(1)) < log_accept_prob) {
          b <- b_transform
          c <- c_dash
          log_like <- logl_new
          list_accept_counts$count_accept4 = list_accept_counts$count_accept4 + 1
        }
      }
    }
    
    #************************************
    #DATA AUGMENTATION 
    #************************************
    if (FLAGS_LIST$DATA_AUG){
      
      if(i == 2){
        print('DATA AUG CHECK')
      }
      #FOR EACH S_T
      for(t in 1:time){
        
        #Copy of data (or update as necessary)
        data_dash = data

        #STOCHASTIC PROPOSAL for s
        if (runif(1) < 0.5) {
          st_dash = data[[2]][t] + 1
        } else {
          st_dash = data[[2]][t] - 1 
        }
        
        #ACCEPTANCE PROBABILITY
        data_dash[[2]][t] = st_dash #s_t = st_dash 
        data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t
        
        #CRITERIA FOR S_T & N_T
        if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
          log_accept_ratio = -Inf 
        } else {
          logl_new = LOG_LIKE_SSI_ORIG(data_dash, a, b, c) 
          log_accept_ratio = logl_new - log_like
        }
        
        #METROPOLIS ACCEPTANCE STEP
        if(log(runif(1)) < log_accept_ratio) {
          data <- data_dash
          log_like <- logl_new
          list_accept_counts$count_accept5 = list_accept_counts$count_accept5 + 1
          vec_accept_da[t] =  vec_accept_da[t] + 1
        }
        
        if(i%%thinning_factor == 0 && i >= burn_in_start){
          non_ss[i_thin, t] = data[[1]][t]
          ss[i_thin, t] = data[[2]][t]
        }
      }
    }
    
    #Loglikelihood Check (Passing - no error)
    if (log_like!=LOG_LIKE_SSI_ORIG(data, a, b, c)){
      print(paste0('log_like: ', log_like, 'LOG_LIKE_CALC(): LOG_LIKE_SSI_ORIG(data, a, b, c)', LOG_LIKE_SSI_ORIG(data, a, b, c)))
    } 
    
    if(i%%thinning_factor == 0 && i >= burn_in_start){
      #POPPULATE MODEL PARAMETERS W/ CURRENT VALUES
      a_vec[i_thin] <- a; b_vec[i_thin] <- b
      c_vec[i_thin] <- c; r0_vec[i_thin] <- a + b*c
      log_like_vec[i_thin] <- log_like
      i_thin = i_thin + 1
    }
    
  }
  
  #BURN IN
  # burn_in_start = burn_in_pc*length(a_vec)
  # a_vec = a_vec[burn_in_start:length(a_vec)]
  # b_vec = b_vec[burn_in_start:length(b_vec)]
  # c_vec = c_vec[burn_in_start:length(c_vec)]
  # log_like_vec = log_like_vec[burn_in_start:length(log_like_vec)]
  # non_ss = non_ss[burn_in_start:nrow(non_ss),]
  # ss = ss[burn_in_start:nrow(ss),]
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1) 
  accept_rate4 = 100*list_accept_counts$count_accept4/(n_mcmc-1)
  accept_rate5 = 100*list_accept_counts$count_accept5/((n_mcmc-1)*time) #i x t
  
  #Acceptance rates 
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5,
                           vec_accept_da = vec_accept_da)
  print(list_accept_rates)
  
  
  return(list(a_vec = a_vec, b_vec = b_vec, c_vec = c_vec, r0_vec = r0_vec,
              log_like_vec = log_like_vec,
              list_accept_rates = list_accept_rates, 
              data = data, mat_count_da = mat_count_da, #13, 14
              non_ss = non_ss, #15
              ss = ss)) #16 
}