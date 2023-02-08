
#************************************************************************
#1. RJMCMC: BASE MODEL VS SSEB MODEL
#************************************************************************
RJMCMC_BASE_SSEB <- function(epidemic_data, n_mcmc,
                            mcmc_inputs = 
                              list(param_starts = list(alpha_start = 0.8, beta_start = 0.1, gamma_start = 10),
                                   alpha_star = 0.4, thinning_factor = 10), 
                            sigma_starts = list(sigma_alpha = 0.3, sigma_beta = 0.03,
                                                sigma_gamma = 3, sigma_bg = 5, sigma_ag = 5),
                            priors_list = list(alpha_prior_exp = c(1, 0), beta_prior_ga = c(10, 2/100),
                                               beta_prior_exp = c(0.1,0),
                                               gamma_prior_ga = c(10, 1), gamma_prior_exp = c(0.1,0)),
                            FLAGS_LIST = list(ADAPTIVE = TRUE, ABG_TRANSFORM = TRUE,
                                              PRIOR = TRUE, BETA_PRIOR_GA = FALSE, GAMMA_PRIOR_GA = FALSE,
                                              THIN = TRUE, alpha_transform = FALSE)) { #?alpha_transform
  
  'Returns MCMC samples of SSEB model parameters (alpha, beta, gamma, r0 = alpha + beta*gamma)
  w/ acceptance rates.
  INCLUDES; ADAPTATION, beta-gamma & alpha-gamma transform'
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  print(paste0('num mcmc iters = ', n_mcmc))
  lambda_vec = get_lambda(epidemic_data)
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_inputs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
  }
  
  #INITIALISE MCMC VECTORS
  alpha_vec <- vector('numeric', mcmc_vec_size); beta_vec <- vector('numeric', mcmc_vec_size)
  gamma_vec <- vector('numeric', mcmc_vec_size); r0_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec <- vector('numeric', mcmc_vec_size);
  
  #INITIALISE MCMC[1]
  alpha_vec[1] <- mcmc_inputs$param_starts$alpha_start; beta_vec[1] <- mcmc_inputs$param_starts$beta_start
  gamma_vec[1] <- mcmc_inputs$param_starts$gamma_start; r0_vec[1] <- alpha_vec[1] + beta_vec[1]*gamma_vec[1]
  log_like_vec[1] <- LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_vec[1], beta_vec[1], gamma_vec[1])
  
  #INITIALISE RUNNING PARAMS
  alpha = alpha_vec[1]; beta = beta_vec[1]; gamma = gamma_vec[1]; log_like = log_like_vec[1]
  
  #SIGMA
  sigma_alpha = sigma_starts$sigma_alpha; sigma_beta = sigma_starts$sigma_beta; 
  sigma_gamma = sigma_starts$sigma_gamma = 3; sigma_bg = sigma_starts$sigma_bg; sigma_ag = sigma_starts$sigma_ag
  
  #SIGMA; INITIALISE FOR ADAPTIVE MCMC
  if (FLAGS_LIST$ADAPTIVE){
    
    #SIGMA
    sigma_alpha_vec <- vector('numeric', mcmc_vec_size); sigma_beta_vec <- vector('numeric', mcmc_vec_size)
    sigma_gamma_vec <- vector('numeric', mcmc_vec_size); sigma_bg_vec <- vector('numeric', mcmc_vec_size)
    #sigma_ag_vec <- vector('numeric', mcmc_vec_size);
    
    #SIGMA; INITIALISE FIRST ELEMENT
    sigma_alpha_vec[1] =  sigma_alpha; sigma_beta_vec[1] = sigma_beta; sigma_gamma_vec[1] =  sigma_gamma
    sigma_bg_vec[1] =  sigma_bg #sigma_ag_vec[1] =  sigma_ag
    
    #SIGMA; List of sigma vectors for each iteration of the MCMC algorithm
    sigma_list = list(sigma_alpha_vec = sigma_alpha_vec, sigma_beta_vec = sigma_beta_vec, sigma_gamma_vec = sigma_gamma_vec,
                      sigma_bg_vec = sigma_bg_vec) #sigma_ag_vec = sigma_ag_vec)
    #Other adaptive parameters
    delta = 1/(mcmc_inputs$alpha_star*(1-mcmc_inputs$alpha_star))
    
  } else {
    
    #SIGMA;sigma vectors for each iteration of the MCMC algorithm
    sigma_list = list(sigma_alpha = sigma_alpha, sigma_beta = sigma_beta,
                      sigma_gamma = sigma_gamma, sigma_bg = sigma_bg)
  }
  
  #INITIALISE: ACCEPTANCE COUNTS
  list_accept_counts = list(count_accept_alpha = 0, count_accept_beta = 0, count_accept_gamma = 0,
                            count_accept_bg = 0, count_accept_m1 = 0,  count_accept_m2 = 0,
                            count_reject_m1 = 0,  count_reject_m2 = 0)

  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {

    if (i%%1000 == 0) {
      
      print(paste0('i = ', i))
    }
    
    #****************************************************** s
    #alpha
    alpha_dash <- alpha + rnorm(1, sd = sigma_alpha)
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }

    #log a
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_dash, beta, gamma)
    log_accept_ratio = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_ratio = log_accept_ratio - alpha_dash + alpha #*Actually this is the Acceptance RATIO. ACCEPTANCE PROB = MIN(1, EXP(ACCPET_PROB))
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      alpha <- alpha_dash
      list_accept_counts$count_accept_alpha = list_accept_counts$count_accept_alpha + 1
      log_like = logl_new
    }
    
    #Sigma (Adaptive)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma_alpha =  sigma_alpha*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #************************************************************************ Only if (b > 0)
    #beta
    beta_dash <- beta + rnorm(1, sd = sigma_beta)
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    #loglikelihood
    logl_new = LOG_LIKE_SSEB(epidemic_data,  lambda_vec, alpha, beta_dash, gamma)
    log_accept_ratio = logl_new - log_like
    
    #Priors
    if (FLAGS_LIST$BETA_PRIOR_GA){
      log_accept_ratio = log_accept_ratio +
        dgamma(beta_dash, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE) -
        dgamma(beta, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - beta_dash + beta
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      beta <- beta_dash
      log_like = logl_new
      list_accept_counts$count_accept_beta = list_accept_counts$count_accept_beta + 1
    }
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma_beta =  sigma_beta*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #************************************************************************
    #gamma
    gamma_dash <- gamma + rnorm(1, sd = sigma_gamma)
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha, beta, gamma_dash)
    log_accept_ratio = logl_new - log_like
    
    #Priors
    if(FLAGS_LIST$GAMMA_PRIOR_GA){
      log_accept_ratio = log_accept_ratio + dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[1], log = TRUE) -
        dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
    } else {
      log_accept_ratio = log_accept_ratio - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
      if (i == 3) print('exp prior on')
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
      gamma <- gamma_dash
      log_like <- logl_new
      list_accept_counts$count_accept_gamma = list_accept_counts$count_accept_gamma + 1
    }
    
    #Sigma (Adpative)
    if (FLAGS_LIST$ADAPTIVE){
      accept_prob = min(1, exp(log_accept_ratio))
      sigma_gamma =  sigma_gamma*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
    }

    #*****************************************************
    #Beta-Gamma TRANSFORM
    if(FLAGS_LIST$ABG_TRANSFORM){
      gamma_dash <- gamma + rnorm(1, sd = sigma_bg)
      
      #Prior > 1 #* TRY WITHOUT REFLECTION
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash
      }
      #New b
      beta_transform = ((alpha + beta*gamma) - alpha)/gamma_dash #beta = (r0 - a)c
      
      if( beta_transform >= 0){ #Only accept values of beta> 0
        
        logl_new = LOG_LIKE_SSEB(epidemic_data,  lambda_vec, alpha, beta_transform, gamma_dash)
        log_accept_ratio = logl_new - log_like
        
        #PRIORS
        #Beta prior
        if (FLAGS_LIST$BETA_PRIOR_GA) {
          tot_beta_prior = dgamma(beta_transform, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE) -
            dgamma(beta, shape = priors_list$beta_prior_ga[1], scale = priors_list$beta_prior_ga[2], log = TRUE)
        } else {
          tot_beta_prior = - beta_transform + beta #exp(1) prior
        }
        
        #gamma prior
        if (FLAGS_LIST$GAMMA_PRIOR_GA) {
          tot_gamma_prior = dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE) -
            dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
        } else {
          tot_gamma_prior = - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
        }
        
        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio + tot_beta_prior + tot_gamma_prior
        
        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          beta <- beta_transform
          gamma <- gamma_dash
          log_like <- logl_new
          list_accept_counts$count_accept_bg = list_accept_counts$count_accept_bg + 1
        }
        
        #Sigma (Adpative)
        if (FLAGS_LIST$ADAPTIVE){
          accept_prob = min(1, exp(log_accept_ratio))
          sigma_bg = sigma_bg*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    }

    #*****************************************************
    #Alpha-Gamma TRANSFORM
    if(FLAGS_LIST$ABG_TRANSFORM){
      
      gamma_dash <- gamma+ rnorm(1, sd = sigma_ag)
      #Prior > 1
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash
      }
      #New alpha
      alpha_transform = (alpha + beta*gamma) - beta*gamma_dash #alpha = (r0 - beta*gamma)
      
      if( alpha_transform >= 0){ #Only accept values of beta> 0
        
        logl_new = LOG_LIKE_SSEB(epidemic_data,  lambda_vec, alpha_transform, beta, gamma_dash)
        log_accept_ratio = logl_new - log_like
        
        #PRIORS
        #gamma prior
        if (FLAGS_LIST$GAMMA_PRIOR_GA) {
          tot_gamma_prior = dgamma(gamma_dash, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE) -
            dgamma(gamma, shape = priors_list$gamma_prior_ga[1], scale = priors_list$gamma_prior_ga[2], log = TRUE)
        } else {
          tot_gamma_prior = - priors_list$gamma_prior_exp[1]*gamma_dash + priors_list$gamma_prior_exp[1]*gamma
        }
        
        #LOG ACCEPT PROB
        log_accept_ratio = log_accept_ratio - alpha_transform + alpha + tot_gamma_prior
        
        #Metropolis Step
        if (!(is.na(log_accept_ratio)) && log(runif(1)) < log_accept_ratio) {
          alpha <- alpha_transform
          gamma <- gamma_dash
          log_like <- logl_new
          list_accept_counts$count_accept_ag = list_accept_counts$count_accept_ag + 1
        }
        
        #Sigma (Adpative)
        if (FLAGS_LIST$ADAPTIVE){
          accept_prob = min(1, exp(log_accept_ratio))
          sigma_ag = sigma_ag*exp(delta/(1+i)*(accept_prob - mcmc_inputs$alpha_star))
        }
      }
    }

    #************************************************************
    #RJMCMC STEP 
    #************************************************************
    
    #************************************************************
    #* M_I *#
    if ((beta > 0) | (gamma > 0)){ #Look to it 
      
      #print('B 0 proposal')
      beta_dash = 0
      gamma_dash = 0
      
      #alpha (alpha only parameter in model)
      if (FLAGS_LIST$alpha_transform) {  #(R0_base) = alpha_sse + beta_sse*gamma_sse (R0_SSE)
        alpha_dash = alpha + beta*gamma #Increase. as alpha_dash is actually the new R_0. Encapsulates 'total R0'          
      } else alpha_dash = alpha #alpha + rnorm(1, sd = sigma_a) - simple stochastic update
      
      #Check alpha positive ==
      if (alpha_dash > 0) { #Automatically satisfied as we've increased alpha. *Remove
        #Acceptance probability (everything cancels)
        logl_new = LOG_LIKE_BASELINE(epidemic_data, alpha_dash)
        
        log_accept_prob = logl_new - log_like - alpha_dash + alpha
        #logl_prev. #Multiply by 100 for example. Increase prior ratio so more likely to accept some M1/spends adequate time in  
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) { #+log(1000)
          beta <- beta_dash
          gamma <- gamma_dash
          alpha <- alpha_dash
          log_like <- logl_new
          list_accept_counts$count_accept_m1 =  list_accept_counts$count_accept_m1 + 1
          
          
        } else  list_accept_counts$count_reject_m1 =  list_accept_counts$count_reject_m1 + 1
      } else  list_accept_counts$count_reject_m1 =  list_accept_counts$count_reject_m1 + 1

    } else { 
      
      #************************************************************
      #* M_II 
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) 
      gamma_dash = rexp(1) + 1 
      
      #alpha
      if (FLAGS_LIST$alpha_transform) { #alpha_sse = ro_base (alpha_base) - beta_sse*gamma_sse
        alpha_dash = alpha - beta_dash*gamma_dash #(alpha_vec[i] - (beta_vec[i]*gamma_vec[i])) - beta_dash*gamma_dash #Preserves alpha, beta, gamma. Will we need the Jacobian?
      } else alpha_dash = alpha #alpha + rnorm(1, sd = sigma_a) # Simple stochastic update
      
      #Check alpha positive==
      if (alpha_dash > 0) {
        
        #Everything cancels 
        logl_new = LOG_LIKE_SSEB(epidemic_data, lambda_vec, alpha_dash, beta_dash, gamma_dash)
        log_accept_prob = logl_new - log_like - alpha_dash + alpha
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta <- beta_dash
          gamma <- gamma_dash
          alpha <- alpha_dash
          log_like <- logl_new
          list_accept_counts$count_accept_m2 =  list_accept_counts$count_accept_m2 + 1
          
        } else list_accept_counts$count_reject_m2 =  list_accept_counts$count_reject_m2 + 1 
      } else list_accept_counts$count_reject_m2 =  list_accept_counts$count_reject_m2 + 1 
    }
    
    #POPULATE VECTORS (ONLY STORE THINNED SAMPLE)
    if (i%%thinning_factor == 0) {
      #print(paste0('i = ', i))
      i_thin = i/thinning_factor
      alpha_vec[i_thin] <- alpha; beta_vec[i_thin] <- beta
      gamma_vec[i_thin] <- gamma; r0_vec[i_thin] <- alpha + beta*gamma
      log_like_vec[i_thin] <- log_like
      sigma_list$sigma_alpha_vec[i_thin] = sigma_alpha; sigma_list$sigma_beta_vec[i_thin] = sigma_beta; 
      sigma_list$sigma_gamma_vec[i_thin] = sigma_gamma
      sigma_list$sigma_bg_vec[i_thin] = sigma_bg; #sigma$sigma_ag_vec[i_thin] = sigma_ag
    }
  }
  
  #Final stats
  accept_rate_a = 100*list_accept_counts$count_accept_alpha/(n_mcmc-1)
  accept_rate_b = 100*list_accept_counts$count_accept_beta/(list_accept_counts$count_accept_beta + list_accept_counts$count_reject_beta)
  accept_rate_g = 100*list_accept_counts$count_accept_gamma/(list_accept_counts$count_accept_gamma + list_accept_counts$count_reject_gamma)
  accept_rate_bg = 100*list_accept_counts$count_accept_bg/(list_accept_counts$count_accept_bg + list_accept_counts$count_reject_bg)
  #RJMCMC Steps 
  accept_rate_m1 = 100*list_accept_counts$count_accept_m1/(list_accept_counts$count_accept_m1 +list_accept_counts$ count_reject_m1) #Check count_accept + count_reject = n_mcmc 
  accept_rate_m2 = 100*list_accept_counts$count_accept_m2/(list_accept_counts$count_accept_m2 + list_accept_counts$count_reject_m2)
  
  #Acceptance rates
  list_accept_rates = list(accept_rate_a = accept_rate_a,
                           accept_rate_b = accept_rate_b, accept_rate_g = accept_rate_g,
                           accept_rate_bg = accept_rate_bg, 
                           accept_rate_m1 = accept_rate_m1, accept_rate_m2 = accept_rate_m2)
  print(list_accept_rates)
  
  #Bayes Factor
  beta_pc0 = (length(which(beta_vec == 0)))/length(beta_vec) 
  bayes_factor = beta_pc0/(1-beta_pc0); bayes_factor = round(bayes_factor, 4)
  
  #OUTPUT
  mcmc_output = list(alpha_vec = alpha_vec, beta_vec = beta_vec, gamma_vec = gamma_vec, r0_vec = r0_vec,
                     log_like_vec = log_like_vec, sigma_list = sigma_list,
                     list_accept_rates = list_accept_rates, beta_pc0 = beta_pc0, bayes_factor = bayes_factor)
  
  #saveRDS(mcmc_output, file = 'mcmc_sse_output_poisson_compound.rds')
  
  return(mcmc_output)
}



