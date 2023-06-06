#MODEL EVIDENCE FOR SSI MODEL
library(extraDistr)
library(MultiRNG)

#FUNCTIONS
ddirmnom(r_dir_samps, 3, alpha.vec, log = FALSE)
rdirmnom(n, size, alpha)


#R SAMP DIRICHLET  MULTINOMIAL
R_MULTINOM_DIR_SS_PROSOSAL <- function(x, mcmc_output, num_is_samps = 1000, beta = 1){
  
  #PARAMS
  N = dim(mcmc_output$ss_matrix)[1] #Sum of the counts of each category, i.e num of 0s + num of 1s (I.e num of mcmc runs)
  print(paste0('N', N))
  matrix_rmultdir_samps = matrix(0, nrow = num_is_samps, ncol = length(x)) #ncol = time
  
  for (t in 1:length(x)){
    
    print(paste0('mcmc_output$ss_matrix[,t]', length(mcmc_output$ss_matrix[,t])))
    print(table(mcmc_output$ss_matrix[,t]))
    
    categories = unique(mcmc_output$ss_matrix[,t])
    alpha_vec = as.vector(table(mcmc_output$ss_matrix[,t])) #table returns counts of each category 
    print(alpha_vec)
    r_dir = draw.dirichlet.multinomial(no.row = 1, #num_is_samps, 
                                          d = length(categories), #length(alpha_vec),
                                          alpha = alpha_vec, 
                                          beta = beta, #scale 
                                          N = N)#Sum of the counts of each category
    
    r_samp_t = rep(categories, times = r_dir)
    r_samp_t = sample(r_samp_t) #shuffle output
    print(r_samp_t)
    matrix_rmultdir_samps[, t] = r_samp_t
    
  }
  
  return(matrix_rmultdir_samps)
}


DENSITY_PROSOSAL_SS_MULTINOM_DIR <- function(theta_samples_proposal_ss, mcmc_output, alphas){
  
  #FUNCTION
  alpha = length of x_t + 1
  ddirmnom(theta_samples_proposal_ss, size, alpha, log = FALSE)
  
}

GET_LOG_PROPOSAL_Q <- function(mcmc_samples, epidemic_data, FLAGS_MODELS,
                               n_is_samples, dof = 3, prob = 0.95) { #*NEED MCMC_OUTPUT
  
  #PARAMETERS REQUIRED 
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*n_is_samples; 
  samp_size_prior = n_is_samples - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR (FROM PARAMETRIC APPROXIMATION)
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(colMeans(mcmc_samples), each = samp_size_proposal) 
  
  theta_samples_proposal_ss = R_MULTINOM_DIR_SS_PROSOSAL(num_is_samps, epidemic_data, mcmc_output) #size num_is_samps x time
  
  theta_samples_proposal = cbind(theta_samples_proposal, theta_samples_proposal_ss)
  
  #PRIOR 
  theta_samples_prior = GET_PRIOR_THETA_SAMPLES(samp_size_prior, n_dim, FLAGS_MODELS)
  #PRIORS ON S?? 
  
  theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  
  #DEFENSE MIXTURE
  matrix_means =  matrix(rep(colMeans(mcmc_samples), each = n_is_samples), ncol = n_dim)
  
  #DENSITY OF PROPOSAL
  log_proposal_density = dmvt(theta_samples - matrix_means,
                              sigma = cov(mcmc_samples), df = dof) #log = TRUE log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  log_proposal_density_ss = DENSITY_PROSOSAL_SS_MULTINOM_DIR(theta_samples_proposal_ss)
  
  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples,
                                            samp_size_prior, n_dim, FLAGS_MODELS)
  
  #PROPOSAL 
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_priors)) #1 x n_is_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_priors = log_priors)
  
  return(imp_samp_comps)  
}

#APPLY
#APPLY
r_samp_t = R_MULTINOM_DIR_SS_PROSOSAL(data_ssir2, mcmc_output)
r_samp_t

wherby 
data_ssir2 = c(2, 0, 1, 0, 5, 3, 4, 2, 1, 4, 5, 4, 3, 3, 6, 4, 4, 8, 7, 12, 13, 15, 15, 12, 24, 26, 26, 41, 32, 38)
