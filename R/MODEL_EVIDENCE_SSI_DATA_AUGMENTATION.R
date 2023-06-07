#MODEL EVIDENCE FOR SSI MODEL
library(SuperSpreadingEpidemicsMCMC)
library(extraDistr)
library(MultiRNG)

#FUNCTIONS
ddirmnom(r_dir_samps, 3, alpha.vec, log = FALSE)
rdirmnom(n, size, alpha)

#OTHER
# DENSITY_PROSOSAL_SS_MULTINOM_DIR <- function(theta_samples_proposal_ss, mcmc_output, alphas){
#   
#   #FUNCTION
#   alpha = length of x_t + 1
#   ddirmnom(theta_samples_proposal_ss, size, alpha, log = FALSE)
#   
# }

#R SAMP DIRICHLET  MULTINOMIAL
PROSOSAL_SS_DIR_MULTINOM <- function(x, mcmc_output, num_is_samps = 1000, beta = 1){
  
  #PARAMS
  N = dim(mcmc_output$ss_matrix)[1] #Sum of the counts of each category, i.e num of 0s + num of 1s (I.e num of mcmc runs)
  matrix_rdirmult_samps = matrix(0, nrow = num_is_samps, ncol = length(x)) #ncol = time
  density_dirmult_samps = c()
    
  for (t in 1:length(x)){
    
    print(paste0('t = ', t))
    
    categories = unique(mcmc_output$ss_matrix[,t])
    alpha_vec = as.vector(table(mcmc_output$ss_matrix[,t])) #table returns counts of each category 
    r_dir_multinom = rdirmnom(n = 1, size = N, alpha = alpha_vec)
    
    # r_dir = draw.dirichlet.multinomial(no.row = 1, #num_is_samps, 
    #                                       d = length(categories), #length(alpha_vec),
    #                                       alpha = alpha_vec, 
    #                                       beta = beta, #scale 
    #                                       N = N)#Sum of the counts of each category
    # 
    # print('r_dir 1')
    # print(r_dir)
    
    r_samp_t = rep(categories, times = r_dir_multinom)
    r_samp_t = sample(r_samp_t) #shuffle output

    matrix_rdirmult_samps[, t] = r_samp_t
    
    density_dirmult_samps[t] = ddirmnom(r_dir_multinom, N, alpha_vec, log = TRUE) #Question should this be applied to r_sampt or r_dir. r_dir matches with alpha
    
  }
  
  return(list(matrix_rdirmult_samps = matrix_rdirmult_samps, density_dirmult_samps = density_dirmult_samps))
}

GET_LOG_PROPOSAL_Q_SS <- function(mcmc_samples, epidemic_data,
                                  num_is_samps, dof = 3, prob = 0.95) { #*NEED MCMC_OUTPUT
  
  #EXTRACT
  dir_multinom_comps = PROSOSAL_SS_DIR_MULTINOM(epidemic_data, mcmc_output, num_is_samps)
  theta_samples_proposal_ss = dir_multinom_comps$matrix_rdirmult_samps   
  density_dirmult_samps = dir_multinom_comps$density_dirmult_samps   
  
  #PARAMS
  lambda_vec = get_lambda(epidemic_data)
  sum_estimate = 0
  
  #SAMPLING SIZE 
  samp_size_proposal = prob*num_is_samps; 
  samp_size_prior = num_is_samps - samp_size_proposal
  prob_prop = prob; prob_prior = 1 - prob_prop
  
  #THETA SAMPLES: PROPOSAL + PRIOR (FROM PARAMETRIC APPROXIMATION)
  theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
    rep(colMeans(mcmc_samples), each = samp_size_proposal)
  
  theta_samples_prior = GET_PRIOR_THETA_SAMPLES(samp_size_prior, n_dim) #FLAGS_MODELS)
  
  #theta_samples_proposal = cbind(theta_samples_proposal, theta_samples_proposal_ss)
  #PRIOR 
  #theta_samples = rbind(theta_samples_proposal, theta_samples_prior)
  
  #DENSITY OF PROPOSAL
  matrix_means =  matrix(rep(colMeans(mcmc_samples), each = n_is_samples), ncol = n_dim)
  log_proposal_density = dmvt(theta_samples - matrix_means,
                              sigma = cov(mcmc_samples), df = dof) #log = TRUE log of the density of multi-variate t distribution (if x = 1,  y= 2, f(x,y) = -4.52) for examples
  
  #PRIOR DENSITIES 
  log_prior_density = GET_LOG_PRIOR_DENSITY(theta_samples,
                                            samp_size_prior, n_dim, FLAGS_MODELS)
  
  #log_proposal_density_ss = DENSITY_PROSOSAL_SS_MULTINOM_DIR(theta_samples_proposal_ss)
  
  #DEFENSE MIXTURE 
  log_q = log(prob_prop*exp(log_proposal_density) + prob_prior*exp(log_priors)) #1 x n_is_samples
  
  imp_samp_comps = list(theta_samples = theta_samples, log_q = log_q, log_priors = log_priors)
  
  return(imp_samp_comps)  
}



#**************************************************************************************
#*
#* 2. GET ESTIMATE OF MODEL EVIDENCE (LOG) (P_HAT)
#*
#*****************************************************************************************
GET_LOG_MODEL_EVIDENCE_SSIB <- function(mcmc_output, epidemic_data, n_samples = 1000
                                        FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                            SSIB = TRUE, SSIR = FALSE)) {   
  
  'Estimate of model evidence for SSEB model using Importance Sampling'
  
  #PARAMS
  vector_estimate_terms = rep(NA, n_samples)
  lambda_vec = get_lambda(epidemic_data); 
  
  #PROPOSAL, PRIOR, THETA SAMPLES 
  mcmc_param_samples = matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
  imp_samp_comps = GET_LOG_PROPOSAL_Q(mcmc_param_samples, epidemic_data, FLAGS_MODELS, n_samples)
  theta_samples = imp_samp_comps$theta_samples
  log_q = imp_samp_comps$log_q; log_prior_density = imp_samp_comps$log_prior_density
  
  #SS MULTINOM DIR
  dir_multinom_comps = PROSOSAL_SS_DIR_MULTINOM(epidemic_data, mcmc_output, num_is_samps)
  theta_samples_proposal_ss = dir_multinom_comps$matrix_rdirmult_samps   
  log_density_dirmult_samps = dir_multinom_comps$density_dirmult_samps   
  log_prior_so = log(1/(1 + epidemic_data[1]))
  
    #GET ESTIMATE
    for (i in 1:n_samples) {
      if (i %% 100 == 0)
        print(i)
      
      loglike = LOG_LIKE_DATA_AUG_SSIB(epidemic_data, theta_samples_proposal_ss[i,], theta_samples[i, 1],
                              theta_samples[i, 2], theta_samples[i, 3]) #theta_samples_proposal_ss
      
      vector_estimate_terms[i] = loglike + log_prior_density[i] + log_prior_so -
        log_q[i] - log_density_dirmult_samps[i]
    }
    
  log_model_ev_est = -log(n_samples) + LOG_SUM_EXP(vector_estimate_terms)
  print(paste0('log_model_ev_est = ', log_model_ev_est))
  
  
  return(log_model_ev_est)
  
}

#APPLY

#MOCK DATA
data_ssir2 = c(2, 0, 1, 0, 5, 3, 4, 2, 1, 4, 5, 4, 3, 3, 6, 4, 4, 8, 7, 12, 13, 15, 15, 12, 24, 26, 26, 41, 32, 38)
n_data = length(data_ssir2)
mcmc_output = list()
ss_matrix = matrix(round(runif(1000, 0, 3)), nrow = 100, ncol = n_data)
mcmc_output$ss_matrix = ss_matrix 

#APPLY
r_samp_t = R_MULTINOM_DIR_SS_PROSOSAL(data_ssir2, mcmc_output)
r_samp_t

help(runif)

a = rdirmnom(1000, 10, c(1:10))
