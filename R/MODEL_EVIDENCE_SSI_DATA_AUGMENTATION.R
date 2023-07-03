#MODEL EVIDENCE FOR SSI MODEL

# - Dirichlet multinomial Importance sampling proposal for SSIB model with Data Augmentation - 
#library(SuperSpreadingEpidemicsMCMC)
#library(mvtnorm)
library(extraDistr)
#library(MultiRNG)

#FUNCTIONS
#ddirmnom(r_dir_samps, 3, alpha.vec, log = FALSE)
#rdirmnom(10, 1, c(2,3))

#R SAMP DIRICHLET  MULTINOMIAL
#' Title
#'
#' @param x 
#' @param mcmc_output 
#' @param num_is_samps 
#' @param beta 
#'
#' @return
#' @export
#'
#' @examples

#**************************************************************************************
#*
#* 1. PROPOSALS FROM DIRICHLET MULTINOMIAL
#*
#*****************************************************************************************
PROSOSAL_SS_DIR_MULTINOM <- function(x, mcmc_output, num_is_samps = 1000, beta = 0.1){ #beta strictly less than 1 
  
  #PARAMS
  N = dim(mcmc_output$ss)[1] #num of mcmc runs     #Sum of the counts of each category, i.e num of 0s + num of 1s (I.e )
  matrix_rdirmult_samps = matrix(0, nrow = num_is_samps, ncol = length(x)) #ncol = time
  density_dirmult_samps = rep(0, num_is_samps) #c()
  print(paste0('beta = ', beta))
    
  for (t in 1:length(x)){
    
    categories = sort(unique(mcmc_output$ss[,t]))
    
    alpha_vec = as.vector(table(mcmc_output$ss[,t])) #table returns counts of each category 
    alpha_vec = alpha_vec*beta
    #print(paste0('alpha_vec = ', alpha_vec))
    
    if (length(alpha_vec)== 1){
      
      matrix_rdirmult_samps[, t] = rep(categories, num_is_samps) #Contributes zero probability 
    } else {
    
      r_dir_multinom = rdirmnom(n = num_is_samps, size = 1, alpha = alpha_vec) #Binary matrix
      
      r_samp_t = r_dir_multinom%*%categories
      
      matrix_rdirmult_samps[, t] = r_samp_t

      density_dirmult_samps = density_dirmult_samps + ddirmnom(x = r_dir_multinom, size = 1, alpha = alpha_vec, log = TRUE) 
    }
    
  }
  
  #print('***********')
  #print(r_dir_multinom)
  #print('***********')
  #print(r_samp_t)
  
  return(list(matrix_rdirmult_samps = matrix_rdirmult_samps, density_dirmult_samps = density_dirmult_samps))
} #beta = 1

#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE 
#*
#***************************************************************************************

LOG_LIKE_DATA_AUG_SSIB <- function(epidemic_data, ss, aX, bX, cX,
                                   shape_gamma = 6, scale_gamma = 1){
  
  #Data
  non_ss = epidemic_data - ss
  
  if(min(non_ss)<0){
    print('WARNING')
    print(non_ss)
    print(ss)
  }
  #non_ss = pmax(non_ss, 0)
  
  #Params
  num_days = length(epidemic_data)
  loglike = 0
  
  #INFECTIOUSNESS  - Difference of two GAMMA distributions. Discretized
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) { 
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS
    lambda_t = sum((non_ss[1:(t-1)] + cX*ss[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    #LOG-LIKELIHOOD
    loglike_t = - lambda_t*(aX + bX) + non_ss[t]*(log(aX) + log(lambda_t)) +
      ss[t]*(log(bX) + log(lambda_t)) - lfactorial(non_ss[t]) - lfactorial(ss[t])
    
    loglike = loglike + loglike_t
    
  }
  
  return(loglike)
}

#' @export
GET_LOG_MODEL_EVIDENCE_SSIB <- function(mcmc_output, epidemic_data, num_is_samps = 1000,
                                        beta = 0.1, FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                            SSIB = TRUE, SSIR = FALSE)) {   
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  num_is_samps = length(mcmc_output$a_vec) #QUESTION 3
  vector_estimate_terms = rep(NA, num_is_samps)
  lambda_vec = get_lambda(epidemic_data) 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  mcmc_param_samples = matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_param_samples, epidemic_data, FLAGS_MODELS, num_is_samps)
  theta_samples = imp_samp_comps$theta_samples
  #print(paste0('theta_samps:', theta_samples))
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #SS MULTINOM DIR
  dir_multinom_comps = PROSOSAL_SS_DIR_MULTINOM(epidemic_data, mcmc_output, num_is_samps, beta = beta)
  theta_samples_proposal_ss = dir_multinom_comps$matrix_rdirmult_samps   
  log_density_dirmult_samps = dir_multinom_comps$density_dirmult_samps   
  
  log_prior_so = log(1/(1 + epidemic_data[1]))
  
    #GET ESTIMATE
    for (i in 1:num_is_samps) {
      
      if (log_prior_density[i] > -Inf){
      
      loglike = LOG_LIKE_DATA_AUG_SSIB(epidemic_data, theta_samples_proposal_ss[i,], theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3]) #theta_samples_proposal_ss
      
      vector_estimate_terms[i] = loglike + log_prior_density[i] + log_prior_so -
        log_q[i] - log_density_dirmult_samps[i]
      
      } else {
        vector_estimate_terms[i] =  -Inf
      }
    }
  
  #print(vector_estimate_terms)
  log_model_ev_est = -log(num_is_samps) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_model_ev_est = ', log_model_ev_est))
  
  
  return(log_model_ev_est)
  
}

#PROPOSAL CHECK 
#dir_multi_nom_comps = PROSOSAL_SS_DIR_MULTINOM(data_ssib, mcmc_output)
#dir_multi_nom_comps


#APPLY
#model_ev_ssib = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_ssib, data_ssib) #mcmc_output

#r_dir_multinom = rdirmnom(n = 1, size = N, alpha = alpha_vec)
# r_dir_multinom = draw.dirichlet.multinomial(no.row = 1, #num_is_samps,
#                                       d = length(categories), #length(alpha_vec),
#                                       alpha = alpha_vec,
#                                       beta = beta, #scale
#                                       N = N)#Sum of the counts of each category