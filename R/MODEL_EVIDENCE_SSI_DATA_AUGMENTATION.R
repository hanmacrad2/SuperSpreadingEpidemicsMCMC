#Model Evidence for SSI model with data augmentation
library(MultiRNG)

#BACKGROUND
#- The Dirichlet distribution is used with the Multinomial when we consider a Bayesian approach to modeling the distribution. In simple words, the probabilities 
# are the parameters of the counts 
# and the priors 
# are the parameters of the probabilities

#SSIB
#SSIB MODEL
#Simulate from a multi-variate multi-dirichelt dist for each timepoint (each column)
#Columns are timepoints, rows are importance samples 
#One importance sample
#Outer loop is over the importance samples 
#One loop over timepoints. Each entry is an interger == number of super-spreaders
#Proposal -> multi-nomial 


#NEED
theta_samples_proposal = rmvt(samp_size_proposal, sigma = cov(mcmc_samples), df = dof) +
  rep(means, each = samp_size_proposal) 

log_proposal_density = dmvt(theta_samples - matrix_means,
                            sigma = cov(mcmc_samples), df = dof) 

#SETUP
ss = matrix(round(runif(1000, 1, 20)), nrow = 100, ncol = 10)

#Generation
#MULTI-NOMIAL DIRICHELT
alpha = round(colSums(non_ss)*0.01)
N = sum(alpha)
beta = 1
no.row =  dim(non_ss)[1]
d =  dim(non_ss)[2]
draw.dirichlet.multinomial(no.row, d, alpha, beta, N) #no.row; (1st element)

r_dir = draw.dirichlet.multinomial(no.row, d, alpha, beta, N) # :D


#Density function (mass function) for the multi-nomial dirichlet model
MULTI_DIRICHLET <- function(data_vector, shrink_factor = 0.1){
  
  data_vector = shrink_factor*data_vector
  term1 = factorial(sum(data_vector)-1)
  term2 = prod(factorial(data_vector-1))
  term3 = prod(data_vector^(data_vector-1))
  
  prob = (term1/term2)*term3
  
  return(prob)
}

#TOY EXAMPLE
data_ssi = c(200, 600, 200)

MULTI_DIRICHLET(data_ssi)

MULTI_DIRICHLET <- function(data_vector, shrink_factor = 0.1) {
  data_vector <- shrink_factor * data_vector
  
  log_term1 <- lgamma(sum(data_vector))
  log_term2 <- sum(lgamma(data_vector))
  log_term3 <- sum((data_vector - 1) * log(data_vector))
  
  log_prob <- log_term1 - log_term2 + log_term3
  prob <- exp(log_prob)
  
  return(prob)
}

# Example usage
data_ssi <- c(200, 600, 200)
prob <- MULTI_DIRICHLET(data_ssi)
print(prob)
