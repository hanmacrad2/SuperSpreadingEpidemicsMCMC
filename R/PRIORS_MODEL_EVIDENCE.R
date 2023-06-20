#MODEL EVIDENCE: FUNCTIONS FOR PRIORS

#DEFINE PRIORS
SET_PRIORS <- function(list_priors = list(priors_sseb = list(exp_prior = c(1,0)),
                                          priors_ssnb = list(pk_prior_nb = c(1,0),
                                                             pk_ga_shape = 0.001, pk_ga_rte = 0.001,
                                                             pr0_unif = c(1.0,4), p_prob_unif = c(0,1)),
                                          priors_ssir = list(pk_exp = c(1,0), pR0_exp = c(1,0)),
                                          priors_ssib = list(exp_prior = c(1,0))),
                       PRIORS_USED = list(SSNB_K_EXP = TRUE, SSNB_K_GAMMA = FALSE)) {
  
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
    pr0_unif = list_priors$priors_ssnb$pr0_unif 
    p_prob_unif =  list_priors$priors_ssnb$p_prob_unif 
    
    if (PRIORS_USED$SSNB_K_EXP){
    
      theta_samples_prior = matrix(c(rexp(samp_size_prior), runif(samp_size_prior,  min = pr0_unif[1],
                                                                  max = pr0_unif[2])), ncol = n_dim)
      
    } else if (PRIORS_USED$SSNB_K_GAMMA){ 
      
      #PRIORS
      pk_ga_shape = list_priors$priors_ssnb$pk_ga_shape
      pk_ga_rte = list_priors$priors_ssnb$pk_ga_rte
      
      theta_samples_prior = matrix(c(rgamma(samp_size_prior, shape = pk_ga_shape, rate = pk_ga_rte),
                                     runif(samp_size_prior,  min = pr0_unif[1], max = pr0_unif[2])),
                                   ncol = n_dim)
    }
    
  } else if (FLAGS_MODELS$SSIR) {
    
    #PRIORS
    pk_exp = list_priors$priors_ssnb$pk_exp
    pR0_exp = list_priors$priors_ssnb$pR0_exp
    
    param_priors = cbind(rexp(samp_size_prior, rate = pk_exp[1]),
                         rexp(samp_size_prior, rate = pR0_exp[1]))
    
    eta_priors_matrix = GET_SAMPLES_ETA_PRIORS(param_priors, epidemic_data, samp_size_prior)
    
    theta_samples_prior = matrix(c(param_priors, eta_priors_matrix), ncol = n_dim)
     
  } else if (FLAGS_MODELS$SSIB) {
    
    #PRIORS
    param_priors = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior),
                                   (1 + rexp(samp_size_prior))), ncol = n_dim) 
    
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
GET_LOG_PRIOR_DENSITY <- function(theta_samples, epidemic_data, samp_size_prior, n_dim, FLAGS_MODELS){
  
  #PRIORS
  list_priors = SET_PRIORS()$list_priors
  PRIORS_USED =  SET_PRIORS()$PRIORS_USED
  
  #1. SSEB
  if(FLAGS_MODELS$SSEB){
    
    log_prior_density = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE) 
    
    #2. SSNB
  } else if (FLAGS_MODELS$SSNB){
    pr0_unif = list_priors$priors_ssnb$pr0_unif 
    p_prob_unif =  list_priors$priors_ssnb$p_prob_unif 
    
    if (PRIORS_USED$SSNB_K_EXP){
      log_prior_density = dexp(theta_samples[,1], log = TRUE) 
      dunif(theta_samples[, 2], min = pr0_unif[1], max = pr0_unif[2], log = TRUE)
    } else if (PRIORS_USED$SSNB_K_GAMMA) {
      log_prior_density = dgamma(theta_samples[,1], shape = list_priors$priors_ssnb$pk_ga_shape, rate = list_priors$priors_ssnb$pk_ga_rte, log = TRUE) +
        dunif(theta_samples[, 2], min =pr0_unif[1], max = pr0_unif[2], log = TRUE) 
    }
    
  } else if (FLAGS_MODELS$SSIR) {
    
    #PRIORS
    log_density_eta_priors = GET_DENSITY_ETA_PRIORS(theta_samples, epidemic_data)
    log_prior_density = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) + log_density_eta_priors
    
  } else if (FLAGS_MODELS$SSIB){
    
    log_prior_density = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE) 
    
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
  
  time_length =  length(epidemic_data) - 1
  eta_priors_matrix = matrix(nrow = samp_size_prior, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size_prior){
    R0X = param_priors[i, 1]; k = mcmc_samples[i, 2]
    eta_priors_matrix[i, ] = rgamma(time_length, shape = epidemic_data*k, scale = R0X*k)
  }
  
  return(eta_priors_matrix)
}

#************************
#DENSITY
#************************
GET_DENSITY_ETA_PRIORS <- function(theta_samples, epidemic_data){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  
  num_etas =  length(epidemic_data) - 1
  samp_size = dim(theta_samples)[1]
  dim_cols = dim(theta_samples)[2] - 1 
  
  eta_samples_matrix = matrix(nrow = samp_size, ncol = num_etas)
  
  #For each mcmc run
  for(i in 1:samp_size){
    R0 = theta_samples[i, 1]; k = mcmc_samples[i, 2]
    
    #print(paste0('dim2: ', length(theta_samples[i, 3:dim_cols])))
    density_samples = dgamma(theta_samples[i, 3:dim_cols], 
                             shape = epidemic_data[1:num_etas]*k,
                             scale = R0*k, log = TRUE)
    
    eta_samples_matrix[i, ] = density_samples #dgamma(theta_samples[i,3:dim_cols], shape = epidemic_data*k, scale = R0X*k, log = TRUE)
  }
  
  return(eta_samples_matrix)
}


