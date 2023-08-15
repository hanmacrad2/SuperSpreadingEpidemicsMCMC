#***********
#PRIORS
#***********

#SSE
GET_PRIORS_SSE <- function() {
  
  PRIORS_SSE = list(R0 = c(1,0), #exp
         k =  c(1,0)) #exp
  
  return(PRIORS_SSE)
}

#SSI
GET_PRIORS_SSI <- function(){
  
  PRIORS_SSI = list(R0 = c(1,0), #exp
                    k =  c(1,0)) #exp
  
  return(PRIORS_SSI)
}

#SSE-B
GET_PRIORS_SSEB <- function() {
  
  PRIORS_SSEB = list(R0 = c(1,0),    #exp dist 
                    alpha =  c(1,2), #beta dist
                    gamma = c(2,1)) #gamma dist
  
  return(PRIORS_SSE)
}

#SSI-B
GET_PRIORS_SSIB <- function(){
  
  PRIORS_SSIB = list(R0 = c(1,0),    #exp dist 
                     a =  c(1,2), #beta dist
                     c = c(2,1)) #gamma dist
  
  return(PRIORS_SSIB)
}

#GET PRIORS USED
GET_PRIORS_USED <- function(){
  
  PRIORS_USED = 
    list(BASELINE = 
         list(R0 = list(EXP = TRUE, GAMMA = FALSE)),
         SSE = 
        list(R0 = list(EXP = TRUE, UNIF = FALSE),
             k =  list(EXP = TRUE, GAMMA = FALSE)),
        SSI = 
          list(R0 = list(EXP = TRUE),
               k =  list(EXP = TRUE)),
        SSEB =
          list(R0 = list(EXP = TRUE),
               alpha = list(EXP = FALSE, BETA = TRUE),
               gamma = list(EXP = FALSE, GAMMA = TRUE))
        SSIB = 
          list(R0 = list(EXP = TRUE),
               a = list(EXP = FALSE, BETA = TRUE),
               c = list(EXP = FALSE, GAMMA = TRUE)))

  return(PRIORS_USED)
}

SET_PRIORS <- function(){
  
  return(list(list_priors = list(
    priors_sse = GET_PRIORS_SSE(),
    priors_ssi = GET_PRIORS_SSI(),
    priors_sseb = GET_PRIORS_SSEB(),
    priors_ssib = GET_PRIORS_SSIB()
  ), PRIORS_USED = GET_PRIORS_USED()))
}

#DEFINE PRIORS
SET_PRIORS <- function(list_priors = list(priors_sseb = list(exp_prior = c(1,0)),
                                          priors_sse = list(pk_prior_nb = c(1,0),
                                                             pk_ga_shape = 0.001, pk_ga_rte = 0.001,
                                                             pr0_unif = c(1.0,4), p_prob_unif = c(0,1)),
                                          priors_ssir = list(pk_exp = c(1,0), pR0_exp = c(1,0)),
                                          priors_ssib = list(a_prior_exp = c(1, 0),
                                                             b_prior_exp = c(1,0),
                                                             c_prior_exp = c(0.1,0),
                                                             a_prior_gamma = c(2, 0.6))),
                       PRIORS_USED = list(BASELINE_EXP = TRUE, BASELINE_GAMMA = FALSE,
                                          SSE_R0_EXP = TRUE, SSE_K_EXP = TRUE, 
                                          SSE_K_GAMMA = FALSE, SSE_R0_UNIF = FALSE,
                                          SSIB_GAMMA_A = FALSE)) {
  
  return(list(list_priors = list_priors, PRIORS_USED = PRIORS_USED))
}

#******************************** 
#* 
#* #1. SAMPLES FROM PRIORS 
#* 
#********************************
GET_PRIOR_THETA_SAMPLES <- function(epidemic_data, samp_size_prior, n_dim, FLAGS_MODELS){
 
  #PRIORS
  list_priors = SET_PRIORS()$list_priors
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  
  #1. SSEB
  if(FLAGS_MODELS$SSEB){
    #PRIORS
    theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = n_dim) 
    
    #2. SSNB
  } else if (FLAGS_MODELS$SSNB){
    n_dim = 2 #CHECK 
    
    if (PRIORS_USED$SSE_R0_EXP && PRIORS_USED$SSE_K_EXP){
    
      theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior)), ncol = n_dim)
    
     } else if (PRIORS_USED$SSE_R0_UNIF){
       
       pr0_unif = list_priors$priors_sse$pr0_unif 
       p_prob_unif =  list_priors$priors_sse$p_prob_unif 
       theta_samples_prior = matrix(c(rexp(samp_size_prior), runif(samp_size_prior,  min = pr0_unif[1],
                                                                   max = pr0_unif[2])), ncol = n_dim)
        
    } else if (PRIORS_USED$SSE_K_GAMMA){ 
      
      #PRIORS
      pk_ga_shape = list_priors$priors_sse$pk_ga_shape
      pk_ga_rte = list_priors$priors_sse$pk_ga_rte
      
      theta_samples_prior = matrix(c(rgamma(samp_size_prior, shape = pk_ga_shape, rate = pk_ga_rte),
                                     runif(samp_size_prior,  min = pr0_unif[1], max = pr0_unif[2])),
                                   ncol = n_dim)
    }
    
  } else if (FLAGS_MODELS$SSIR) {
    
    #PRIORS
    #browser()
    pR0_exp = list_priors$priors_ssir$pR0_exp
    pk_exp = list_priors$priors_ssir$pk_exp
    
    prior_R0 =  rexp(samp_size_prior,  rate = pR0_exp[1]) 
    #problem_R0 = which(prior_R0 < 0.09)
    #prior_R0[problem_R0] = prior_R0[problem_R0] + runif(length(prior_R0[problem_R0]), 0.1, 0.5)
    
    prior_k =  rexp(samp_size_prior,  rate = pk_exp[1]) 
    #problem_k = which(prior_k < 0.09)
    #prior_k[problem_k] = prior_k[problem_k] + runif(length(prior_k[problem_k]), 0.1, 0.5)
    
    param_priors = cbind(prior_R0, prior_k)
    
    # param_priors = cbind(rexp(samp_size_prior, rate = pR0_exp[1]),
    #                      rexp(samp_size_prior, rate = pk_exp[1]))
    
    #ETAS
    eta_priors_matrix = GET_SAMPLES_ETA_PRIORS(param_priors, epidemic_data, samp_size_prior)
    theta_samples_prior = cbind(param_priors, eta_priors_matrix)
    #theta_samples_prior = matrix(c(param_priors, eta_priors_matrix), ncol = n_dim): 
    print(paste0('dim theta_samples_prior', dim(theta_samples_prior)))
     
  } else if (FLAGS_MODELS$SSIB) {
    
    #PRIORS
    if (PRIORS_USED$SSIB_GAMMA_A) {
      
      gamma_shape = list_priors$priors_ssib$a_prior_gamma[1]
      gamma_scale = list_priors$priors_ssib$a_prior_gamma[2]
      param_priors = matrix(c(rgamma(samp_size_prior, shape = gamma_shape, scale = gamma_scale), rexp(samp_size_prior),
                              (1 + rexp(samp_size_prior))), ncol = n_dim) 
    } else {
      param_priors = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior),
                              (1 + rexp(samp_size_prior, rate = list_priors$priors_ssib$c_prior_exp[1]))), ncol = n_dim) 
    }
    
    theta_samples_prior = param_priors
    
    #data_aug_priors_matrix = GET_SAMPLES_DATA_AUG_PRIORS(param_priors, epidemic_data, samp_size_prior)
    #ncol = n_dim + dim_data_aug
    #theta_samples_prior = matrix(c(param_priors, eta_priors_matrix), ncol = n_dim)
    
  }
  
  return(theta_samples_prior)
}

#********************************
#* 
#* 2. GET DENSITY (LOG) OF PRIORS 
#*
#********************************
GET_LOG_PRIOR_DENSITY <- function(theta_samples, epidemic_data,
                                  samp_size_prior, n_dim, FLAGS_MODELS){
  
  #browser()
  #PRIORS
  list_priors = SET_PRIORS()$list_priors
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  
  #1. SSEB
  if(FLAGS_MODELS$SSEB){
    
    log_prior_density = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE) 
    
    #2. SSNB
  } else if (FLAGS_MODELS$SSNB){
    
    if (PRIORS_USED$SSE_K_EXP && PRIORS_USED$SSE_R0_EXP){
      
      log_prior_density = dexp(theta_samples[,1], log = TRUE)  + 
        dexp(theta_samples[, 2], log = TRUE) 
      
    } else if (PRIORS_USED$SSE_K_GAMMA) {
      pr0_unif = list_priors$priors_sse$pr0_unif 
      p_prob_unif =  list_priors$priors_sse$p_prob_unif 
      log_prior_density = dgamma(theta_samples[,1], shape = list_priors$priors_sse$pk_ga_shape, rate = list_priors$priors_sse$pk_ga_rte, log = TRUE) +
        dunif(theta_samples[, 2], min =pr0_unif[1], max = pr0_unif[2], log = TRUE) 
    }
    
  } else if (FLAGS_MODELS$SSIR) {
    
    #PRIORS
    pR0_exp = list_priors$priors_ssir$pR0_exp
    pk_exp = list_priors$priors_ssir$pk_exp
    log_density_eta_priors = GET_DENSITY_ETA_PRIORS(theta_samples, epidemic_data)
    log_prior_density = dexp(theta_samples[,1], rate = pR0_exp[1], log = TRUE) +
      dexp(theta_samples[,2], rate = pk_exp[1], log = TRUE) + log_density_eta_priors
    
  } else if (FLAGS_MODELS$SSIB){
    
    #PRIORS
    if (PRIORS_USED$SSIB_GAMMA_A) {
      
      gamma_shape = list_priors$priors_ssib$a_prior_gamma[1]
      gamma_scale = list_priors$priors_ssib$a_prior_gamma[2]
      log_prior_density = dgamma(theta_samples[,1], shape = gamma_shape, scale = gamma_scale, log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
        dexp((theta_samples[,3] - 1), log = TRUE) 
    } else {
      log_prior_density = dexp(theta_samples[,1], log = TRUE) +
        dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), rate = list_priors$priors_ssib$c_prior_exp, log = TRUE) 
    }
    
  }
  
  return(log_prior_density)
}


#************************
#*
#FUNCTIONS FOR DATA AUGMENTATION MODELS (SSIR #1)
#* 
#************************

GET_SAMPLES_ETA_PRIORS <- function(param_priors, epidemic_data, samp_size_prior){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
 
  time_length =  length(epidemic_data) - 1 #ETA LENGTH TIME - 1
  wh_nonzero = which(epidemic_data[1:(length(epidemic_data)-1)]!= 0)
  eta_priors_matrix = matrix(0, nrow = samp_size_prior, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size_prior){
    
    R0 = param_priors[i, 1]; k = param_priors[i, 2]
    eta_prior_vec = c()
    
    for (j in 1:time_length){
      eta_prior = rgamma(1, shape = epidemic_data[j]*k, scale = R0/k)
      # 
      # if (eta_prior < 0.005) { #((eta_prior > 0) && (eta_prior < 0.005))
      #   #browser()
      #   eta_prior = eta_prior +  runif(1, 0.01, 0.2)
      # }
      eta_prior_vec[j] = eta_prior
    }
    #eta_prior = rgamma(time_length, shape = epidemic_data[1:time_length]*k, scale = R0/k)
    
    eta_priors_matrix[i, ] = eta_prior_vec 
  }
  
  return(eta_priors_matrix)
}

#*ORIGAL
V0_GET_SAMPLES_ETA_PRIORS <- function(param_priors, epidemic_data, samp_size_prior){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  #browser()
  time_length =  length(epidemic_data) - 1 #ETA LENGTH TIME - 1
  eta_priors_matrix = matrix(0, nrow = samp_size_prior, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size_prior){ 
    R0 = param_priors[i, 1]; k = param_priors[i, 2]
    eta_prior = rgamma(time_length, shape = epidemic_data[1:time_length]*k, scale = R0/k)
    eta_priors_matrix[i, ] = eta_prior # rgamma(time_length, shape = epidemic_data[1:time_length]*k, scale = R0/k)
    
  }
  
  return(eta_priors_matrix)
}

#************************
#DENSITY
#************************
GET_DENSITY_ETA_PRIORS <- function(theta_samples, epidemic_data){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  #browser()
  samp_size = dim(theta_samples)[1]
  dim_cols = dim(theta_samples)[2] 
  log_density_eta = rep(NA, times = samp_size) 
  
  #For each mcmc run
  for(i in 1:samp_size){
    
    R0 = theta_samples[i, 1]; k = theta_samples[i, 2]
    
    #ADDED 14/07/23
    if (any(theta_samples[i, ] < 0)) { #WARNING CAUSED WITHOUT THIS 
      log_density_eta[i] = -Inf #rep(-Inf, time = num_etas)
    } else {
      wh_nonzero = which(epidemic_data[1:(length(epidemic_data)-1)]!= 0)
      log_density_eta[i] = sum(dgamma(theta_samples[i, 2+wh_nonzero, drop = FALSE], 
                               shape = epidemic_data[wh_nonzero, drop = FALSE]*k, 
                               scale = R0/k, log = TRUE))
    }
    
    #log_density_eta[i] = density_samples 
  }
  
  return(log_density_eta)
}

