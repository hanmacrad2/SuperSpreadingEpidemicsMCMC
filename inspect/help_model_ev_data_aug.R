
#A COLUMN FOR THE PRIOR FOR EACH ETA_T (EVERY POINT IN TIME)
GET_SAMPLES_ETA_PRIORS <- function(param_priors, epidemic_data, samp_size_prior){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  
  time_length =  length(epidemic_data)
  eta_priors_matrix = matrix(nrow = samp_size_prior, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size_prior){
    R0X = param_priors[i, 1]; k = mcmc_samples[i, 2]
    eta_priors_matrix[i, ] = rgamma(time_length, shape = epidemic_data*k, scale = R0X*k)
  }

  return(eta_priors_matrix)
}

GET_DENSITY_ETA_PRIORS <- function(theta_samples, epidemic_data){
  
  'Get priors for all etas in the SSI model'
  #Question: Is it correct to use the current priors r0x & k
  
  time_length =  length(epidemic_data)
  samp_size = dim(theta_samples)[1]
  eta_samples_matrix = matrix(nrow = samp_size, ncol = time_length)
  
  #For each mcmc run
  for(i in 1:samp_size){
    R0X = theta_samples[i, 1]; k = mcmc_samples[i, 2]
    eta_samples_matrix[i, ] = dgamma(time_length, shape = epidemic_data*k, scale = R0X*k, log = TRUE)
  }
  
  return(eta_samples_matrix)
}

#RUN
mcmc_samples = mcmc_samples[,1:11]
epidemic_data = epidemic_data[1:9]
n_dim = 11
means = colMeans(mcmc_samples)
#PART 1 
theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
  rep(means, each = samp_size_proposal) 

#PRIOR
param_priors = cbind(rexp(samp_size_prior),
                     rexp(samp_size_prior))

eta_priors_matrix = GET_PRIORS_ETA(param_priors, epidemic_data, samp_size_prior)

theta_samples_prior = matrix(c(param_priors, eta_priors_matrix), ncol = n_dim)

#Proposal samples
theta_samples = rbind(theta_samples_proposal, theta_samples_prior)

#print(paste0('mean_mcmc = ', means))
print(paste0('mean proposal = ', colMeans(theta_samples_proposal)))
print(paste0('mean prior = ', colMeans(theta_samples_prior)))

#DEFENSE MIXTURE
n_samples = 120
log_proposal = dmvt(theta_samples - matrix(rep(means, each = n_samples), ncol = n_dim),
                    sigma = cov(mcmc_samples), df = dof, log = TRUE)

log_density_eta_priors = GET_DENSITY_ETA_PRIORS(theta_samples, epidemic_data)
log_priors = dexp(theta_samples[,1], log = TRUE) + dexp(theta_samples[,2], log = TRUE) + log_density_eta_priors














#RUN 2
params_priors = cbind(rexp(samp_size_prior),
                      rexp(samp_size_prior))

eta_priors_matrix = GET_PRIORS_ETA(params_priors, data_wait_08_21_sub1, samp_size_prior)

eta_priors_matrix = GET_DENSITY_ETA_PRIORS(params_priors, data_wait_08_21_sub1)

prior_samples = matrix(c(params_priors, eta_priors_matrix), ncol = 12)
       
eta_priors_matrix = GET_PRIORS_ETA(params_priors, epidemic_data, samp_size_prior)

theta_samples_prior = matrix(c(params_priors, eta_priors_matrix))


#COMBINE
n_dim = dim(out_ssir)[2]
eta_priors_matrix = GET_PRIORS_ETA(mcmc_samples, data_wait_08_21_sub1, samp_size_prior)
theta_samples_prior = matrix(c(rexp(samp_size_prior, rate = priors_ssir$pk_exp[1]),
                               rexp(samp_size_prior, rate = priors_ssir$pR0_exp[1]), eta_priors_matrix), ncol = n_dim)



#PRIOR ON ETA = GAMMAS 
eta_prob = dgamma(eta[t-1], shape = x[t-1]*k, scale = R0/k, log = TRUE)
priors_eta = GET_PRIORS_ETA(mcmc_samples, epidemic_data)

#FOR EACH MCMC SAMPLE
n_mcmc = dim(mcmc_samples)[1]
time =  dim(mcmc_samples)[1] - 2

#LOOP TO GET EACH PRIOR (CAN YOU VECTORISE)?
for (i in 1:n_mcmc) { 
  
  #eta_prior = 
  #OUTER Is for the number of samples (vectorised below)
  #theta_samples_prior = matrix(c(rexp(samp_size_prior, rate = priors_ssir$pk_exp[1]),
                                 #rexp(samp_size_prior, rate = priors_ssir$pR0_exp[1]), priors_eta), ncol = n_dim)
  
  #t or t-1??: #NEED TO THINK ABOUT MODEL
  for (t in 1:time) { 
    eta_prior = rgamma(eta[t-1], shape = x[t-1]*k, scale = R0/k)

  }
  
  
}

eta_priors_matrix = matrix(nrow = n_mcmc, ncol = length(epidemic_data))
for(i in 1:n_mcmc){
  
  R0X = mcmc_samples[i, 1]; k = mcmc_samples[i, 2]
  eta_priors_matrix[i, ] = rgamma(length(epidemic_data), shape = epidemic_data*k, scale = R0X*k)
}

for (t in 1:length(x)){ #For the etas (So only simulating one value from gamma at a time) *And compare to vectorised results
  
}

#TRIAL
samp_size_proposal = 100; dof = 3; samp_size_prior = 20
n_dim = 3
mcmc_output = mcmc_sseb
mcmc_samples =  matrix(c(mcmc_output$alpha_vec, mcmc_output$beta_vec, mcmc_output$gamma_vec), ncol = 3)

theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
  rep(colMeans(mcmc_samples), each = samp_size_proposal) 

theta_samples_prior = matrix(c(rexp(samp_size_prior), rexp(samp_size_prior), (1 + rexp(samp_size_prior))), ncol = n_dim)

theta_samples_prior = matrix(c(rgamma(samp_size_prior, shape = 0.001, rate = 0.001),
                               runif(samp_size_prior,  min = 0.5, max = 4), ncol = n_dim))

data_wait_08_21_sub1




#CHECK
a = cbind(c(1,2,3),c(4,5,6))

b =  cbind(c(7,8,9),c(10,11,12))

c = matrix(c(a,b), ncol = 4)


#DATA
shape_gamma = 6; scale_gamma = 1; 
num_days = 10
pgamma(c(1:10), shape = shape_gamma, scale = scale_gamma)

prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
  pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)


#SSIR NAN VALUES
eta_prob = NA; eta_prob2 = -Inf 
if (!is.na(eta_prob) && !is.infinite(eta_prob)){
  print('yes')
} 
