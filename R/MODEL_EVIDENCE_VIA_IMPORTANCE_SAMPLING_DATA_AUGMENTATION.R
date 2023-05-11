#'Estimate of Model Evidence via Importance sampling as in:
#'Touloupou, Panayiota, et al.
#"Efficient model comparison techniques for models requiring large scale data augmentation." (2018): 437-459.'

#LIBRARIES
#library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
library(compositions)

#***************************************
#*
#1. GET IMPORTANCE SAMPLING PROPOSAL (LOG)
#*
#***************************************

#MULTI_DIMENSIONAL PROPOSAL
GET_LOG_PROPOSAL_DATA_AUG <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                                         n_samples, dof = 3, prob = 0.95, 
                                         priors_ssnb = list(pk_ga_shape = 0.001, pk_ga_rte = 0.001, pr0_unif = c(1.0,4),
                                                            p_prob_unif = c(0,1), ),
                                      priors_ssir = list(pk_exp = c(1,0), pR0_exp = c(1,0)),
                                         FLAG_PRIORS = list(SSEB_EXP_PRIOR = TRUE)) {               
  
  #PARAMETERS REQUIRED
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_samples; samp_size_prior = n_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR
  means = colMeans(mcmc_samples)
  #theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) + means 
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(means, each = samp_size_proposal) 
  
  #PRIORS
  #print(paste0('FLAGS_MODELS$SSNB', FLAGS_MODELS))
  if(FLAGS_MODELS$SSEB & FLAG_PRIORS$SSEB_EXP_PRIOR){
    n_dim = dim(mcmc_samples)[2] #3 #PUT IN DICTIONARY
    theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = n_dim) 
    
  } else if (FLAGS_MODELS$SSNB){
    n_dim =  dim(mcmc_samples)[2] #2 
    theta_samples_prior = matrix(c(rgamma(samp_size_prior, shape = priors_ssnb$pk_ga_shape, rate = priors_ssnb$pk_ga_rte),
                                   runif(samp_size_prior,  min = priors_ssnb$pr0_unif[1], max = priors_ssnb$pr0_unif[2])), ncol = n_dim)
  } else if (FLAGS_MODELS$SSIR) {
    n_dim = dim(out_ssir)[2]
    priors_eta = GET_PRIORS_ETA(mcmc_samples, epidemic_data)
    theta_samples_prior = matrix(c(rexp(samp_size_prior, rate = priors_ssir$pk_exp[1]),
                                   rexp(samp_size_prior, rate = priors_ssir$pR0_exp[1]), priors_eta), ncol = n_dim)
    
  }
  
  #Proposal samples
  #print(paste0('dim = theta_samples_proposal', dim(theta_samples_proposal)))
  #print(paste0('dim = theta_samples_prior', dim(theta_samples_prior)))
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #print(paste0('mean_mcmc = ', means))
  #print(paste0('mean proposal = ', colMeans(theta_samples_proposal)))
  #print(paste0('mean prior = ', colMeans(theta_samples_prior)))
  
  #DEFENSE MIXTURE
  log_proposal = dmvt(theta_samples - matrix(rep(means, each = n_samples), ncol = n_dim),
                      sigma = cov(mcmc_samples), df = dof, log = TRUE) #log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #PRIORS
  if(FLAGS_MODELS$SSEB  & FLAG_PRIORS$SSEB_EXP_PRIOR){
    log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) +
      dexp((theta_samples[,3] - 1), log = TRUE) 
  } else if (FLAGS_MODELS$SSNB){
    log_priors = dgamma(theta_samples[,1], shape = priors_ssnb$pk_ga_shape, rate = priors_ssnb$pk_ga_rte, log = TRUE) +
      dunif(theta_samples[, 2], min = priors_ssnb$pr0_unif[1], max = priors_ssnb$pr0_unif[2], log = TRUE)
  }
  
  log_q = log(prob_prop*exp(log_proposal) + prob_prior*exp(log_priors)) #1 x n_samples
  #print(paste0('log_q ', log_q))
  
  #max_el = pmax(log(prob_prop) + log_proposal, log(prob_prior) + log_priors)
  #log_q2 = max_el + log(exp(log(prob_prop) + log_proposal - max_el) + exp(log(prob_prior) + log_priors - max_el))
  #log_q_s = LOG_SUM_EXP(log_q2) #LOG SUM EXP OF TWO COMPONENTS - See if the same
  #print(paste0('log_q_s2 ', log_q_s))
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_priors = log_priors)
  
  return(imp_samp_comps)  
}


#******************************************************************
#* 2b. ESTIMATE OF MODEL EVIDENCE (PHAT) FOR SSEB MODEL  
#*******************************************************************
GET_LOG_P_HAT_DATA_AUG <- function(mcmc_samples, epidemic_data,
                          FLAGS_MODELS, n_samples = 1000,
                          priors_ssnb = list(pk_ga_shape = 0.001, pk_ga_rte = 0.001, pr0_unif = c(1.0,4),
                                             p_prob_unif = c(0,1)),
                          FLAG_PRIORS = list(SSEB_EXP_PRIOR = TRUE)) {    
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(NA, n_samples)
  lambda_vec = get_lambda(epidemic_data); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  imp_samp_comps = GET_LOG_PROPOSAL_Q_MULTI_DIM(mcmc_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_priors = imp_samp_comps$log_priors
  
  #SSEB MODEL 
  if (FLAGS_MODELS$SSEB & FLAG_PRIORS$SSEB_EXP_PRIOR) {
    
    #GET ESTIMATE
    for (i in 1:n_samples) {
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_SSEB(epidemic_data, lambda_vec, theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3])
      
      vector_estimate_terms[i] = loglike + log_priors[i] - log_q[i]
    }
    
    #SSNB MODEL
  } else if (FLAGS_MODELS$SSNB) {
    
    for (i in 1:n_samples) {
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_SSNB(epidemic_data, lambda_vec, theta_samples[i,]) 
      
      vector_estimate_terms[i] = loglike + log_priors[i] - log_q[i]
    }
  }
  
  log_p_hat = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_p_hat = ', log_p_hat))
  
  log_p_hat2 = -log(n_samples) + log(sum(exp(vector_estimate_terms)))
  print(paste0('log_p_hat2 = ', log_p_hat2))
  
  return(log_p_hat)
}
