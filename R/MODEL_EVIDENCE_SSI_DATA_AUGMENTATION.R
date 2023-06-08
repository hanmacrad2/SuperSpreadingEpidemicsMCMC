#MODEL EVIDENCE FOR SSI MODEL

# - Dirichlet multinomial Importance sampling proposal for SSIB model with Data Augmentation - 
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)
library(extraDistr)
library(MultiRNG)

#FUNCTIONS
#ddirmnom(r_dir_samps, 3, alpha.vec, log = FALSE)
#rdirmnom(n, size, alpha)

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
PROSOSAL_SS_DIR_MULTINOM <- function(x, mcmc_output, num_is_samps = 1000, beta = 1){
  
  #PARAMS
  N = dim(mcmc_output$ss)[1] #Sum of the counts of each category, i.e num of 0s + num of 1s (I.e num of mcmc runs)
  #matrix_rdirmult_samps = matrix(0, nrow = num_is_samps, ncol = length(x)) #ncol = time
  matrix_rdirmult_samps = matrix(0, nrow = N, ncol = length(x))
  density_dirmult_samps = c()
    
  for (t in 1:length(x)){
    
    categories = unique(mcmc_output$ss[,t])
    alpha_vec = as.vector(table(mcmc_output$ss[,t])) #table returns counts of each category 
    
    r_dir_multinom = draw.dirichlet.multinomial(no.row = 1, #num_is_samps,
                                          d = length(categories), #length(alpha_vec),
                                          alpha = alpha_vec,
                                          beta = beta, #scale
                                          N = N)#Sum of the counts of each category
    
    #r_dir_multinom = rdirmnom(n = 1, size = N, alpha = alpha_vec)
    
    r_samp_t = rep(categories, times = r_dir_multinom)
    r_samp_t = sample(r_samp_t) #shuffle output
    
    matrix_rdirmult_samps[, t] = r_samp_t
    
    #QUESTION 1: Should this be applied to r_samp_t or r_dir_multinom. r_dir_multinom matches with dim of alpha
    density_dirmult_samps[t] = ddirmnom(r_dir_multinom, N, alpha_vec, log = TRUE) 
    
  }
  
  return(list(matrix_rdirmult_samps = matrix_rdirmult_samps, density_dirmult_samps = density_dirmult_samps))
}

#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE 
#*
#*****************************************************************************************
#*
#' @export
GET_LOG_MODEL_EVIDENCE_SSIB <- function(mcmc_output, epidemic_data, num_is_samps = 1000,
                                        FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                            SSIB = TRUE, SSIR = FALSE)) {   
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  num_is_samps = length(mcmc_output$a_vec) #QUESTION 3
  vector_estimate_terms = rep(NA, num_is_samps)
  lambda_vec = get_lambda(epidemic_data) 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  mcmc_param_samples = matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
  print(paste0('dim of mcmc_samps', dim(mcmc_param_samples)))
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_param_samples, epidemic_data, FLAGS_MODELS, num_is_samps)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #SS MULTINOM DIR
  dir_multinom_comps = PROSOSAL_SS_DIR_MULTINOM(epidemic_data, mcmc_output, num_is_samps)
  theta_samples_proposal_ss = dir_multinom_comps$matrix_rdirmult_samps   
  log_density_dirmult_samps = dir_multinom_comps$density_dirmult_samps   
  log_prior_so = log(1/(1 + epidemic_data[1]))
  
  #print(paste0('log_prior_density ', log_prior_density))
  #print(paste0('theta_samples ', theta_samples))
  #print(paste0('theta_samples_proposal_ss ', theta_samples_proposal_ss))
  
    #GET ESTIMATE
    for (i in 1:num_is_samps) {
      
      loglike = LOG_LIKE_DATA_AUG_SSIB(epidemic_data, theta_samples_proposal_ss[i,], theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3]) #theta_samples_proposal_ss
      
      vector_estimate_terms[i] = loglike + log_prior_density[i] + log_prior_so -
        log_q[i] - log_density_dirmult_samps[i]
    }
  
  print(vector_estimate_terms)
  log_model_ev_est = -log(num_is_samps) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_model_ev_est = ', log_model_ev_est))
  
  
  return(log_model_ev_est)
  
}

#PROPOSAL CHECK 
#dir_multi_nom_comps = PROSOSAL_SS_DIR_MULTINOM(data_ssib, mcmc_output)
#dir_multi_nom_comps


#APPLY
model_ev_ssib = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_output, data_ssib)

