#MCMC SSIB OFFSPRING
#********************************************************
#*
#1. SSIB MCMC JOINT NO DATA AUGMENTATION
#* 
#********************************************************
LOGLIKE_SSIB_OFFSPRING <- function(params, data) {
  
  r0 <- params[1]
  a <- params[2]; b <- params[3]
  p = (1 - a)/(1 - a + a*b)
  log_likelihood <- 0
  
  for (z in data) {
    
    # Calculate the two Poisson components
    lambda1 <- a*r0 + (1-a)*r0/b
    lambda2 <- a*b*r0 + (1-a)*r0
    
    # Calculate the log-probabilities for the mixture components
    log_prob1 <- dpois(z, lambda = lambda1, log = TRUE)
    log_prob2 <- dpois(z, lambda = lambda2, log = TRUE)
    
    # Calculate the log of the weighted sum of probabilities
    log_sum_prob <- LOG_SUM_EXP(c(log((1-p)*exp(log_prob1)), log(p*exp(log_prob2))))
    #print(z)
    #print(log_sum_prob)
    
    if(!is.infinite(log_sum_prob)){
      log_likelihood <- log_likelihood + log_sum_prob
    } 
    print(paste0('log_likelihood: ', log_likelihood))
  }
  
  return(log_likelihood) # Return the negative log-likelihood for minimization
}



#' @export 
MCMC_INFER_SSIB_OFFSPRING <- function(offspring_data, n_mcmc, PRIORS_USED = GET_PRIORS_USED(),
                               param_starts = c(1.5, 0.5, 10),
                               mcmc_inputs = list(dim = 3, target_acceptance_rate = 0.4,
                                                  v0 = 100,  thinning_factor = 10,
                                                  burn_in_pc = 0.2)){    
  
  #INITIALIASE DATA
  #non_ss = data$non_ss
  #ss = data$ss
  #epidemic_data = non_ss + ss
  
  #INITIALIASE PARAMS
  #r0_start = GET_R0_INITIAL_MCMC(epidemic_data)
  #param_starts[1] = r0_start
  #time = length(epidemic_data) 
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
  
  #LOG LIKELIHOOD
  log_like_vec <- vector('numeric', mcmc_vec_size)
  log_like_vec[1] <- LOGLIKE_SSIB_OFFSPRING(ssib_params, offspring_data) #, FLAG_NEGBIN_PARAMATERISATION)
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
    
    if(i == 2) print(paste0('i = ', i))
    
    if(i%%1000 == 0) print(paste0('i = ', i))
    
    #PROPOSAL
    ssib_params_dash = c(ssib_params + mvrnorm(1, mu = rep(0, mcmc_inputs$dim), Sigma = scaling*c_star*sigma_i)) 
    
    #POSTIVE ONLY
    if (min(ssib_params_dash - vec_min) >= 0 && ssib_params_dash[2] < 1){ 
      
      #LOG LIKELIHOOD
      logl_new = LOGLIKE_SSIB_OFFSPRING(ssib_params_dash, offspring_data) 
      
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
    
  } #END FOR LOOP
  
  #Final stats
  accept_rate = 100*count_accept/(n_mcmc-1)
  #Return a, acceptance rate
  return(list(ssib_params_matrix = ssib_params_matrix,
              log_like_vec = log_like_vec, scaling_vec = scaling_vec, 
              accept_rate = accept_rate))
} 


# #MCMC SSIB
# n_mcmc = 100000
# freq_sec_cases = c(199, 57, 18, 6, 5, 1, 2,0,0,0,1,1) 
# mcmc_ssib_offspring = MCMC_INFER_SSIB_OFFSPRING(freq_sec_cases, n_mcmc)
# 
# plot.ts(mcmc_ssib_offspring$ssib_params_matrix)
# 
# #BEST FIT
# r0 = mean(mcmc_ssib_offspring$ssib_params_matrix[,1])
# a = mean(mcmc_ssib_offspring$ssib_params_matrix[,2])
# b = mean(mcmc_ssib_offspring$ssib_params_matrix[,3])
# 
# params_ssib = c(r0, a, b)
# 
# #SSIB FIT
# num_offspring <- length(freq_sec_cases)
# x <- 0:(num_offspring - 1)
# ssib_fit <- GET_OFFSPRING_SSIB(x, params_ssib)
# 
# 
# #FIT 1
# #params_ssib1
# #ssib_fit
# #> params_ssib
# #[1] 1.61587271 0.01341351 8.94911074
# 
# #SAVE
# RESULTS_FOLDER = "~/Github/computing/REAL_DATA/4_OFFSPRING_DISTS/OFFSPRING_HK/RESULTS/"
# filename = 'params_ssib_hk1.rds'
# saveRDS(params_ssib, paste0(RESULTS_FOLDER, filename))
# filename = 'ssib_fit_hk1.rds'
# saveRDS(ssib_fit, paste0(RESULTS_FOLDER, filename))
# filename = 'mcmc_ssib_fit_hk1.rds'
# saveRDS(mcmc_ssib_offspring, paste0(RESULTS_FOLDER, filename))
# 
# file_name = 'mcmc_ssib_offspring_sp.rds'
# saveRDS(mcmc_ssib_offspring_sp, paste0(RESULTS_FOLDER, filename))
# 
# filename = 'ssib_fit_sp.rds'
# saveRDS(ssib_fit_sp, paste0(RESULTS_FOLDER, filename))
# 
# #LOAD MCMC
# mcmc_ssib_offspring = readRDS(paste0(RESULTS_FOLDER, filename))
#   
# #FIT 2
# params_ssib2 = c(r0, a, b)
# 
# #SSIB FIT
# num_offspring <- length(freq_sec_cases)
# x <- 0:(num_offspring - 1)
# ssib_fit2 <- GET_OFFSPRING_SSIB(x, params_ssib)
# 
# 
# #*******************************
# #SINGAPORE DATA
# n_mcmc = 100000
# sec_cases_sp <- c(160, 20, 6, 5, 0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,
#                   1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
# mcmc_ssib_offspring = MCMC_INFER_SSIB_OFFSPRING(freq_sec_cases, n_mcmc)
# 
# plot.ts(mcmc_ssib_offspring$ssib_params_matrix)
# 
# #BEST FIT
# r0 = mean(mcmc_ssib_offspring$ssib_params_matrix[,1])
# a = mean(mcmc_ssib_offspring$ssib_params_matrix[,2])
# b = mean(mcmc_ssib_offspring$ssib_params_matrix[,3])
