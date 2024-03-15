#************************************
#*
#* GET PRIOR IMPORTANCE SAMPLES: 
#* 
#* IMPORTANCE SAMPLES FROM PRIOR DISTRIBUTIONS FOR MIXTURE DISTRIBUTION
#*
#*************************************

#SSEB MODEL
GET_PRIOR_SAMPS_SSEB <- function(samp_size_prior, n_dim = 3){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSEB()
  
  #ALPHA: HOW TO INCLUDE r0
  if (PRIORS_USED$SSEB$alpha$BETA) {
    
    shape1 = list_priors$alpha[1]
    shape2 = list_priors$alpha[2]
    samps_prior_alpha = rbeta(samp_size_prior, shape1, shape2)
  }
  
  #r0
  if (PRIORS_USED$SSEB$r0$EXP) {
    samps_prior_r0 = rexp(samp_size_prior, rate = list_priors$r0[1])
  }
  
  #GAMMA
  if (PRIORS_USED$SSEB$beta$GAMMA) {
    shape = list_priors$beta[1]
    scale = list_priors$beta[2]
    samps_prior_beta = 1 + rgamma(samp_size_prior, shape = shape, scale = scale)
    
  }
  
  #EXP
  # if(PRIORS_USED$SSEB$alpha$EXP){
  #   theta_samples_prior = matrix(c(rexp(samp_size_prior),
  #                                  rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = n_dim) 
  #   
  # }
  
  #THETA_SAMPLES_PRIOR
  theta_samples_prior_sseb = matrix(c(samps_prior_alpha, samps_prior_r0, samps_prior_beta), ncol = n_dim)
  
  return(theta_samples_prior_sseb)  
}

#SSIB MODEL
GET_PRIOR_SAMPS_SSIB <- function(samp_size_prior, n_dim = 3){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSIB()
  
  #r0
  if (PRIORS_USED$SSIB$r0$EXP) {
    samps_prior_r0 = rexp(samp_size_prior, rate = list_priors$r0[1])
  }
  
  #a 
  if (PRIORS_USED$SSIB$a$BETA) {
    shape1 = list_priors$a[1]
    shape2 = list_priors$a[2]
    samps_prior_a = rbeta(samp_size_prior, shape1, shape2)
  }
  
  #b
  if (PRIORS_USED$SSIB$b$GAMMA) {
    shape = list_priors$b[1]
    scale = list_priors$b[2]
    samps_prior_b = 1 + rgamma(samp_size_prior, shape = shape, scale = scale)
    
  }
  
  #THETA_SAMPLES_PRIOR
  theta_samples_prior_ssib = matrix(c(samps_prior_r0, samps_prior_a, samps_prior_b),
                                    ncol = n_dim)
  
  return(theta_samples_prior_ssib)  
}

#SSE MODEL
GET_PRIOR_SAMPS_SSE <- function(samp_size_prior, n_dim = 2){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSE()
  
  #r0
  if (PRIORS_USED$SSE$r0$EXP) {
    samps_prior_r0 = rexp(samp_size_prior, rate = list_priors$r0[1])
  }
  
  #k
  if (PRIORS_USED$SSE$k$EXP) {
    samps_prior_k = rexp(samp_size_prior, rate = list_priors$k[1])
  }
  
  #THETA_SAMPLES_PRIOR
  theta_samples_prior_sse = matrix(c(samps_prior_r0, samps_prior_k), ncol = n_dim)
  
  return(theta_samples_prior_sse)  
}

#SSI MODEL
GET_PRIOR_SAMPS_SSI <- function(samp_size_prior, n_dim = 2){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSI()
  
  #r0
  if (PRIORS_USED$SSI$r0$EXP) {
    samps_prior_r0 = rexp(samp_size_prior, rate = list_priors$r0[1])
  }
  
  #k
  if (PRIORS_USED$SSI$k$EXP) {
    samps_prior_k = rexp(samp_size_prior, rate = list_priors$k[1])
  }
  
  #THETA_SAMPLES_PRIOR
  theta_samples_prior_ssi = cbind(samps_prior_r0, samps_prior_k)
  
  return(theta_samples_prior_ssi)  
}


#****************************************************
#*
#* 3. GET PRIOR DENSITES (LOG) OF IMP SAMPS 
#*
#*****************************************************

#SSEB MODEL
GET_DENSITY_PRIOR_SSEB <- function(theta_samples){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSEB()
  
  #ALPHA 
  if (PRIORS_USED$SSEB$alpha$BETA) {
    shape1 = list_priors$alpha[1]
    shape2 = list_priors$alpha[2]
    log_prior_density_alpha = dbeta(theta_samples[,1], shape1, shape2, log = TRUE)
  }
  
  #r0
  if (PRIORS_USED$SSEB$r0$EXP) {
    log_prior_density_r0 = dexp(theta_samples[,2], rate = list_priors$r0[1], log = TRUE)
  }
  
  #GAMMA
  if (PRIORS_USED$SSEB$beta$GAMMA) {
    shape = list_priors$beta[1]
    scale = list_priors$beta[2]
    log_prior_density_beta = dgamma(theta_samples[,3]-1, shape = shape, scale = scale, log = TRUE)
    
  }
  
  #EXP
  # if(PRIORS_USED$SSEB$alpha$EXP){
  #   log_prior_density = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
  #     dexp((theta_samples[,3] - 1), log = TRUE) 
  # }
  
  #THETA_SAMPLES_PRIOR
  log_prior_density = log_prior_density_alpha + log_prior_density_r0 + log_prior_density_beta
  
  return(log_prior_density)  
}


#SSEB MODEL
GET_DENSITY_PRIOR_SSIB <- function(theta_samples){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSIB()
  
  
  #r0
  if (PRIORS_USED$SSIB$r0$EXP) {
    log_prior_density_r0 = dexp(theta_samples[,1], rate = list_priors$r0[1], log = TRUE)
  }
  
  #a
  if (PRIORS_USED$SSIB$a$BETA) {
    shape1 = list_priors$a[1]
    shape2 = list_priors$a[2]
    log_prior_density_a = dbeta(theta_samples[,2], shape1, shape2, log = TRUE)
  }
  
  #b
  if (PRIORS_USED$SSIB$b$GAMMA) {
    shape = list_priors$b[1]
    scale = list_priors$b[2]
    log_prior_density_b = dgamma(theta_samples[,3]-1, shape = shape, scale = scale, log = TRUE)
    
  }
  
  #THETA_SAMPLES_PRIOR
  log_prior_density = log_prior_density_r0 + log_prior_density_a + log_prior_density_b
  
  return(log_prior_density)  
}

#SSE MODEL
GET_DENSITY_PRIOR_SSE <- function(theta_samples){
  
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSE()
  
  #r0
  if (PRIORS_USED$SSE$r0$EXP) {
    log_prior_density_r0 = dexp(theta_samples[,1], rate = list_priors$r0[1], log = TRUE)
  }
  
  #k
  if (PRIORS_USED$SSE$k$EXP) {
    log_prior_density_k = dexp(theta_samples[,2], rate = list_priors$k[1], log = TRUE)
  }
  
  #THETA_SAMPLES_PRIOR
  log_prior_density =  log_prior_density_r0 + log_prior_density_k
  
  return(log_prior_density)  
}


#SSE MODEL
GET_DENSITY_PRIOR_SSI <- function(theta_samples){
  
  #browser()
  #List priors
  PRIORS_USED = GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSI()
  
  #r0
  if (PRIORS_USED$SSI$r0$EXP) {
    log_prior_density_r0 = dexp(theta_samples[,1], rate = list_priors$r0[1], log = TRUE)
  }
  
  #k
  if (PRIORS_USED$SSI$k$EXP) {
    log_prior_density_k = dexp(theta_samples[,2], rate = list_priors$k[1], log = TRUE)
  }
  
  #THETA_SAMPLES_PRIOR
  log_prior_density =  log_prior_density_r0 + log_prior_density_k
  
  return(log_prior_density)  
}


#******************************** 
#* 
#* #2b. SAMPLES FROM PRIORS 
#* 
#********************************
GET_PRIOR_IMPORTANCE_SAMPLES <- function(epidemic_data, samp_size_prior, FLAGS_MODELS){
  
  #1. SSEB
  if(FLAGS_MODELS$SSEB){
    #PRIORS
    theta_samples_prior = GET_PRIOR_SAMPS_SSEB(samp_size_prior)
    
    #2. SSE
  } else if (FLAGS_MODELS$SSE){
    
    #PRIORS 
    theta_samples_prior = GET_PRIOR_SAMPS_SSE(samp_size_prior)
    
  } else if (FLAGS_MODELS$SSI) {
    
    #PRIORS
    param_priors = GET_PRIOR_SAMPS_SSI(samp_size_prior)
    eta_priors_matrix = GET_SAMPLES_ETA_PRIORS(param_priors, epidemic_data, samp_size_prior)
    theta_samples_prior = cbind(param_priors, eta_priors_matrix)
    
  } else if (FLAGS_MODELS$SSIB) {
    
    #PRIORS
    theta_samples_prior = GET_PRIOR_SAMPS_SSIB(samp_size_prior)
    
  }
  
  return(theta_samples_prior)
}

#********************************
#* 
#* 3b. GET PRIOR DENSITES (LOG) OF IMP SAMPS  
#*
#********************************
GET_LOG_PRIOR_DENSITY <- function(theta_samples, epidemic_data, FLAGS_MODELS){
  
  #1. SSEB
  if(FLAGS_MODELS$SSEB){
  
    log_prior_density = GET_DENSITY_PRIOR_SSEB(theta_samples)
    
    #2. SSE
  } else if (FLAGS_MODELS$SSE){
    
    log_prior_density = GET_DENSITY_PRIOR_SSE(theta_samples)
    
  } else if (FLAGS_MODELS$SSI) {
    #browser()
    log_density_params = GET_DENSITY_PRIOR_SSI(theta_samples)
    log_density_eta_priors = GET_DENSITY_ETA_PRIORS(theta_samples, epidemic_data)
    log_prior_density = log_density_params + log_density_eta_priors
    
  } else if (FLAGS_MODELS$SSIB) {
    
    log_prior_density = GET_DENSITY_PRIOR_SSIB(theta_samples)
  }
  
  return(log_prior_density)
}


#************************************************************************************
#*
#   ETA PRIORS
#
# - FUNCTIONS FOR DATA AUGMENTATION MODEL (SSI)
#* 
#************************************************************************************

GET_SAMPLES_ETA_PRIORS <- function(param_priors, epidemic_data, samp_size_prior){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  
  time_length =  length(epidemic_data) - 1 #ETA LENGTH TIME - 1
  wh_nonzero = which(epidemic_data[1:(length(epidemic_data)-1)]!= 0)
  eta_priors_matrix = matrix(0, nrow = samp_size_prior, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size_prior){
    
    r0 = param_priors[i, 1]; k = param_priors[i, 2]
    eta_prior_vec = c()
    
    for (j in 1:time_length){
      eta_prior = rgamma(1, shape = epidemic_data[j]*k, scale = r0/k)
      eta_prior_vec[j] = eta_prior
    }
    #eta_prior = rgamma(time_length, shape = epidemic_data[1:time_length]*k, scale = r0/k)
    
    eta_priors_matrix[i, ] = eta_prior_vec 
  }
  
  return(eta_priors_matrix)
}

GET_DENSITY_ETA_PRIORS <- function(theta_samples, epidemic_data){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  #browser()
  samp_size = dim(theta_samples)[1]
  dim_cols = dim(theta_samples)[2] 
  log_density_eta = rep(NA, times = samp_size) 
  
  #For each mcmc run
  for(i in 1:samp_size){
    
    r0 = theta_samples[i, 1]; k = theta_samples[i, 2]
    
    if (any(theta_samples[i, ] < 0)) { #WARNING CAUSED WITHOUT THIS 
      log_density_eta[i] = -Inf #rep(-Inf, time = num_etas)
    } else {
      wh_nonzero = which(epidemic_data[1:(length(epidemic_data)-1)]!= 0)
      log_density_eta[i] = sum(dgamma(theta_samples[i, 2+wh_nonzero, drop = FALSE], 
                                      shape = epidemic_data[wh_nonzero, drop = FALSE]*k, 
                                      scale = r0/k, log = TRUE))
    }
    
    #log_density_eta[i] = density_samples 
  }
  
  return(log_density_eta)
}
