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
  mcmc_param_samples = matrix(c(mcmc_output$alpha_vec, mcmc_output$r0_vec, mcmc_output$c_vec), ncol = 3)
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_param_samples, epidemic_data, FLAGS_MODELS, num_is_samps)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #SS MULTINOM DIR
  dir_multinom_comps = PROSOSAL_SS_DIR_MULTINOM(epidemic_data, mcmc_output, num_is_samps, beta = beta)
  theta_samples_proposal_ss = dir_multinom_comps$matrix_rdirmult_samps   
  log_density_dirmult_samps = dir_multinom_comps$density_dirmult_samps   
  
  log_prior_so = log(1/(1 + epidemic_data[1])) #What is this? s_0?
  
  #GET ESTIMATE
  for (i in 1:num_is_samps) {
    
    if (log_prior_density[i] > -Inf){
      
      loglike = LOG_LIKE_DATA_AUG_SSIB(epidemic_data, theta_samples_proposal_ss[i,], theta_samples[i, 1],
                                       theta_samples[i, 2], theta_samples[i, 3]) #theta_samples_proposal_ss
      
    } else {
      loglike = 0 #vector_estimate_terms[i] =  0 #-Inf
    }
    vector_estimate_terms[i] = loglike + log_prior_density[i] + log_prior_so -
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


#VECTORISED
PROSOSAL_SS_DIR_MULTINOM <- function(x, mcmc_output, num_is_samps = 10000,
                                        beta = 1000, prior_dir = 0.8){ #beta strictly less than 1 #prior dir: try 0 too
  
  #PARAMS
  N = dim(mcmc_output$ss)[1] #num of mcmc runs     #Sum of the counts of each category, i.e num of 0s + num of 1s (I.e )
  
  matrix_rdirmult_samps = matrix(NA, nrow = num_is_samps, ncol = length(x)) #ncol = time
  
  density_dirmult_samps = rep(0, num_is_samps) #c()
  
  #print(paste0('beta = ', beta))
  
  for (t in 1:length(x)){
    
    if (t == 1){
      last_categories = c(NA) 
    } else {
      #num of SS proposed at previous time point
      last_categories = sort(unique(matrix_rdirmult_samps[,t-1])) #from all is at last timepoint. find categories #categories 
    }
    
    for (j in 1:length(last_categories)) {
      
      if (t == 1) {
        wh_mcmc = 1:N #Take all mcmc samples for first time point  
        wh_is = 1:num_is_samps #importance samples we're going to sample. Variable depending on num of categories in prev. iteration
      } else { #Matches previous time point
        
        wh_mcmc = which(mcmc_output$ss[,t-1] == last_categories[j]) #matrix_rdirmult_samps[last_categories, t-1]) #otherwise find mcmc iterations that match the importance sample at the previous time point
        wh_is = which(matrix_rdirmult_samps[,t-1] == last_categories[j]) # last_categories[j] == the category
      }
      
      categories = sort(unique(mcmc_output$ss[wh_mcmc,t])) #Take the subset where mcmc == import. samp
      
      #Match at previous time point and see where it goes. Conditional dist of st|st-1. Simulate from the part 
      alpha_vec = as.vector(table(mcmc_output$ss[wh_mcmc,t])) #table returns counts of each category 
      #Normalise
      alpha_vec = (alpha_vec*beta)/sum(alpha_vec) #+ rep(prior_dir/length(categories), length(categories))  #sum(alpha) = effective sample size of the prior
      #beta = effective sample size #Jim Burger: prior_dir = 0.8
      
      if (length(alpha_vec)== 1){ #in case matrix == one dimension 
        
        matrix_rdirmult_samps[wh_is, t] = categories #rep(categories, num_is_samps) #Contributes zero log probability. As only one category so prob = 1, log(prob) = 0.  
        #density_dirmult_samps[wh_is] = density_dirmult_samps[wh_is] + 0: Does nothing. 
        
      } else {
        
        r_dir_multinom = rdirmnom(n=length(wh_is), size = 1, alpha = alpha_vec) #Binary matrix
        
        r_samp_t = r_dir_multinom%*%categories
        
        matrix_rdirmult_samps[wh_is, t] = r_samp_t
        
        density_dirmult_samps[wh_is] = density_dirmult_samps[wh_is] + ddirmnom(x = r_dir_multinom, size = 1, alpha = alpha_vec, log = TRUE) 
      }
      
    }
    
  }
  
  return(list(matrix_rdirmult_samps = matrix_rdirmult_samps, density_dirmult_samps = density_dirmult_samps))
} 


#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE 
#*
#***************************************************************************************

LOG_LIKE_DATA_AUG_SSIB <- function(epidemic_data, ss, alpha, r0, c,
                                   shape_gamma = 6, scale_gamma = 1){
  
  #Params
  b = (r0*(1 - alpha))/c #r0 = a_prop*r0 + b*c
  a = alpha*r0 
  
  #Data
  non_ss = epidemic_data - ss
  
  # if(min(non_ss)<0){
  #   print('WARNING')
  #   print(non_ss)
  #   print(ss)
  # }
  
  #Params
  num_days = length(epidemic_data)
  loglike = 0
  
  #INFECTIOUSNESS  - Difference of two GAMMA distributions. Discretized
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) { 
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS
    lambda_t = sum((non_ss[1:(t-1)] + c*ss[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    #LOG-LIKELIHOOD
    loglike_t = dpois(non_ss[t], a*lambda_t, log = TRUE) +
      dpois(ss[t], b*lambda_t, log = TRUE)
    
    # loglike_t = - lambda_t*(a + b) + non_ss[t]*(log(a) + log(lambda_t)) +
    #   ss[t]*(log(b) + log(lambda_t)) - lfactorial(non_ss[t]) - lfactorial(ss[t])
    
    loglike = loglike + loglike_t
    
  }
  
  return(loglike)
}


#PROPOSAL CHECK 
#dir_multi_nom_comps = PROSOSAL_SS_DIR_MULTINOM(data_ssib, mcmc_output)


#APPLY
#model_ev_ssib = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_ssib, data_ssib) 

#r_dir_multinom = rdirmnom(n = 1, size = N, alpha = alpha_vec)
# r_dir_multinom = draw.dirichlet.multinomial(no.row = 1, #num_is_samps,
#                                       d = length(categories), #length(alpha_vec),
#                                       alpha = alpha_vec,
#                                       beta = beta, #scale
#                                       N = N)#Sum of the counts of each category