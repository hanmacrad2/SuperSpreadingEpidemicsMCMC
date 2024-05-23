#MODEL EVIDENCE FOR SSI MODEL
#Dirichlet multinomial Importance sampling proposal for SSIB model with Data Augmentation - 

#FUNCTIONS
#ddirmnom(r_dir_samps, 3, alpha.vec, log = FALSE)
#rdirmnom(10, 1, c(2,3))


# MODEL EVIDENCE

#' @export
GET_LOG_MODEL_EVIDENCE_SSIB <- function(mcmc_output, epidemic_data, 
                                        num_is_samps = 10000, beta = 1000, 
                                        FLAGS_MODELS = list(BASE = FALSE, SSE = FALSE, SSI = FALSE,
                                                            SSEB = FALSE, SSIB = TRUE)) {   
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  #PARAMS
  vector_estimate_terms = rep(NA, num_is_samps)
  lambda_vec = get_lambda(epidemic_data) 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  mcmc_param_samples = mcmc_output$ssib_params_matrix 
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_param_samples, epidemic_data, FLAGS_MODELS, num_is_samps)
  
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #SS MULTINOM DIR
  dir_multinom_comps = PROSOSAL_SS_DIR_MULTINOM(epidemic_data, mcmc_output, beta = beta)
  theta_samples_proposal_ss = dir_multinom_comps$matrix_rdirmult_samps   
  log_density_dirmult_samps = dir_multinom_comps$density_dirmult_samps   
  
  #log_prior_so = log(1/(1 + epidemic_data[1])) #What is this? s_0?
  
  #GET ESTIMATE
  for (i in 1:num_is_samps) {
    
    if (log_prior_density[i] > -Inf){
      
      loglike = LOG_LIKE_DATA_AUG_SSIB(epidemic_data, theta_samples_proposal_ss[i,], theta_samples[i, ])
                                       #theta_samples[i, 2], theta_samples[i, 3]) #theta_samples_proposal_ss
      
    } else {
      loglike = 0 #vector_estimate_terms[i] =  0 #-Inf
    }
    vector_estimate_terms[i] = loglike + log_prior_density[i] - #+ log_prior_so 
      log_q[i] - log_density_dirmult_samps[i]
  }
  
  log_model_ev_est = -log(num_is_samps) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_model_ev_est = ', log_model_ev_est))
  
  return(log_model_ev_est)
}


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
#' 
#' 

PROSOSAL_SS_DIR_MULTINOM <- function(epidemic_data, mcmc_output, num_is_samps = 10000,
                                     beta = 1000, prior_dir = 0.8) {
  # PARAMS
  N <- dim(mcmc_output$ss_mcmc)[1] # num of mcmc runs
  num_columns <- dim(mcmc_output$ss_mcmc)[2] # number of columns in ss_mcmc matrix
  
  matrix_rdirmult_samps <- matrix(NA, nrow = num_is_samps, ncol = length(epidemic_data)) # ncol = time
  density_dirmult_samps <- rep(0, num_is_samps) # c()
  
  for (t in 1:length(epidemic_data)) {
    if (t == 1) {
      last_categories <- c(NA)
    } else {
      last_categories <- sort(unique(matrix_rdirmult_samps[, t-1]))
    }
    
    for (j in 1:length(last_categories)) {
      if (t == 1) {
        wh_mcmc <- 1:N # Take all mcmc samples for first time point
        wh_is <- 1:num_is_samps # importance samples we're going to sample
      } else {
        if (t-1 > num_columns) {
          #print("t-1 exceeds the number of columns in mcmc_output$ss_mcmc, skipping this iteration")
          next
        }
        wh_mcmc <- which(mcmc_output$ss_mcmc[, t-1] == last_categories[j])
        wh_is <- which(matrix_rdirmult_samps[, t-1] == last_categories[j])
      }
      
      if (length(wh_mcmc) == 0) {
        #print("No matches found in MCMC")
        next
      } else {
        #print(paste("Matches found in MCMC:", length(wh_mcmc)))
      }
      
      if (t > num_columns) {
        print(paste("t exceeds the number of columns in mcmc_output$ss_mcmc, skipping iteration at t =", t))
        next
      }
      
      categories <- sort(unique(mcmc_output$ss_mcmc[wh_mcmc, t]))
      
      alpha_vec <- as.vector(table(mcmc_output$ss_mcmc[wh_mcmc, t]))
      alpha_vec <- alpha_vec * max(beta / sum(alpha_vec), 1) + rep(prior_dir / length(categories), length(categories))
      
      if (length(alpha_vec) == 1) {
        matrix_rdirmult_samps[wh_is, t] <- categories
      } else {
        r_dir_multinom <- rdirmnom(n = length(wh_is), size = 1, alpha = alpha_vec)
        r_samp_t <- r_dir_multinom %*% categories
        matrix_rdirmult_samps[wh_is, t] <- r_samp_t
        density_dirmult_samps[wh_is] <- density_dirmult_samps[wh_is] + ddirmnom(x = r_dir_multinom, size = 1, alpha = alpha_vec, log = TRUE)
      }
    }
  }
  
  return(list(matrix_rdirmult_samps = matrix_rdirmult_samps, density_dirmult_samps = density_dirmult_samps))
}


#VECTORISED
# PROSOSAL_SS_DIR_MULTINOM <- function(x, mcmc_output, num_is_samps = 10000, # beta = 1000
#                                      beta = 1000, prior_dir = 0.8){ #beta strictly less than 1 #prior dir: try 0 too
#   
#   #PARAMS
#   N = dim(mcmc_output$ss_mcmc)[1] #num of mcmc runs     #Sum of the counts of each category, i.e num of 0s + num of 1s (I.e )
#   
#   matrix_rdirmult_samps = matrix(NA, nrow = num_is_samps, ncol = length(x)) #ncol = time
#   
#   density_dirmult_samps = rep(0, num_is_samps) #c()
#   
#   #print(paste0('beta = ', beta))
#   
#   for (t in 1:length(x)){
#     
#     if (t == 1){
#       last_categories = c(NA) 
#     } else {
#       #num of SS proposed at previous time point
#       last_categories = sort(unique(matrix_rdirmult_samps[,t-1])) #from all is at last timepoint. find categories #categories 
#     }
#     
#     for (j in 1:length(last_categories)) {
#       
#       if (t == 1) {
#         wh_mcmc = 1:N #Take all mcmc samples for first time point  
#         wh_is = 1:num_is_samps #importance samples we're going to sample. Variable depending on num of categories in prev. iteration
#       } else { #Matches previous time point
#         
#         wh_mcmc = which(mcmc_output$ss_mcmc[,t-1] == last_categories[j]) #matrix_rdirmult_samps[last_categories, t-1]) #otherwise find mcmc iterations that match the importance sample at the previous time point
#         wh_is = which(matrix_rdirmult_samps[,t-1] == last_categories[j]) # last_categories[j] == the category
#       }
#       
#       if (length(wh_mcmc) == 0) {
#         # Handle the case where there are no matches
#         print("No matches found")
#         browser()
#       } else {
#         # Process the matches
#         print(wh_mcmc)
#       }
#     
#       categories = sort(unique(mcmc_output$ss_mcmc[wh_mcmc, t])) #Take the subset where mcmc == import. samp
#       
#       #Match at previous time point and see where it goes. Conditional dist of st|st-1. Simulate from the part 
#       alpha_vec = as.vector(table(mcmc_output$ss_mcmc[wh_mcmc,t])) #table returns counts of each category 
#       #Normalise
#       #Increase over-dispersion (reduce the tail on model evidence distribution). Outliers greatly increasing variance
#       alpha_vec = alpha_vec*max(beta/sum(alpha_vec), 1) + rep(prior_dir/length(categories), length(categories))
#       #alpha_vec = (alpha_vec*beta)/sum(alpha_vec) #+ rep(prior_dir/length(categories), length(categories))  #sum(alpha) = effective sample size of the prior
#       
#       #beta = effective sample size #Jim Burger: prior_dir = 0.8
#       
#       if (length(alpha_vec)== 1){ #in case matrix == one dimension 
#         
#         matrix_rdirmult_samps[wh_is, t] = categories #rep(categories, num_is_samps) #Contributes zero log probability. As only one category so prob = 1, log(prob) = 0.  
#         #density_dirmult_samps[wh_is] = density_dirmult_samps[wh_is] + 0: Does nothing. 
#         
#       } else {
#         
#         r_dir_multinom = rdirmnom(n=length(wh_is), size = 1, alpha = alpha_vec) #Binary matrix
#         
#         r_samp_t = r_dir_multinom%*%categories
#         
#         matrix_rdirmult_samps[wh_is, t] = r_samp_t
#         
#         density_dirmult_samps[wh_is] = density_dirmult_samps[wh_is] + ddirmnom(x = r_dir_multinom, size = 1, alpha = alpha_vec, log = TRUE) 
#       }
#     }
#     
#   }
#   
#   return(list(matrix_rdirmult_samps = matrix_rdirmult_samps, density_dirmult_samps = density_dirmult_samps))
# } 


#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE 
#*
#***************************************************************************************

LOG_LIKE_DATA_AUG_SSIB <- function(epidemic_data, ss, ssib_params,
                                   shape_gamma = 6, scale_gamma = 1){
  
  #PARAMS
  r0 = ssib_params[1];  a = ssib_params[2]; b = ssib_params[3]
  c = (r0*(1 - a))/b #r0 = a_prop*r0 + b*c
  
  #Data
  non_ss = epidemic_data - ss
  
  #Params
  num_days = length(epidemic_data)
  loglike = 0
  
  #INFECTIOUSNESS  - Difference of two GAMMA distributions. Discretized
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) { 
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS
    lambda_t = sum((non_ss[1:(t-1)] + b*ss[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    #LOG-LIKELIHOOD
    loglike_t = dpois(non_ss[t], a*lambda_t, log = TRUE) +
      dpois(ss[t], c*lambda_t, log = TRUE)
    
    loglike = loglike + loglike_t
    
  }
  
  return(loglike)
}

