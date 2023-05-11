
#A COLUMN FOR THE PRIOR FOR EACH ETA_T (EVERY POINT IN TIME)
GET_PRIORS_ETA <- function (){
  
}

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



for (t in 1:length(x)){ #For the etas (So only simulating one value from gamma at a time) *And compare to vectorised results
  
}