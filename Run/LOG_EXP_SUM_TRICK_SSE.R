#LOG_EXP_SUM TRICK

#1. LOG LIKELIHOOD
LOG_LIKE_SSE_POISSON <- function(x, lambda_vec, alphaX, betaX, gammaX) {
    
  #Params
    num_days = length(x)
    loglike = 0
    #prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
    
    for (t in 2:num_days) {
      #lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)])); inner_sum_xt = 0
      term1_alpha = alphaX*lambda_vec[t]
      
      if ((x[t] == 0) | is.na(x[t])) {
        loglike = loglike - term1_alpha + LSE_ZT(0, lambda_vec[t], betaX, gammaX)
        
      } else {
        inner_sum_vec <- vector('numeric', x[t])
        for (yt in 0:x[t]) {
          #Sum for all values of yt
          
          #Log likelihood
          zt = x[t] - yt
          inner_sum_vec[yt + 1] = -term1_alpha + yt*log(term1_alpha) - lfactorial(yt)
          + LSE_ZT(zt, lambda_vec[t], betaX, gammaX)
        }
        
        #Calculate max element in inner vector for all yts
        lx_max = max(inner_sum_vec)
        
        lse = lx_max + log(sum(exp(inner_sum_vec - lx_max)))
        #print(paste0('1 lse = ', lse))
        if(!is.na(lse)) {
          loglike = loglike + lse 
        }
      }
    }
    
    #print(paste0('loglike', loglike))
    return(loglike)
  }

#LSE GAMMA
LSE_ZT <- function(zt, lambda_t, betaX, gammaX, max_nt = 5) {
  #SETUP
  inner_sum_vec <- vector('numeric', max_nt)
  term1 = betaX*lambda_t
  
  #NT = 0
  inner_sum_vec[1] = -term1
  for (nt in 1:max_nt) {
    term_add = -term1 + nt*log(term1) -gammaX*nt + zt*log(gammaX*nt) -lfactorial(nt)
    
    if(is.na(term_add)){
      print('LSE_ZT is na')
      print(paste0('nt:', nt))
      print(paste0('zt:', nt))
    }
    inner_sum_vec[nt] = term_add
  }
  
  #LSE
  lx_max = max(inner_sum_vec)
  lse = lx_max + log(sum(exp(inner_sum_vec - lx_max)))
  #print(paste0('2 lse = ', lse))

  return(lse)
}

#1. LOG LIKELIHOOD
LOG_LIKE_SSE_POISSON <-
  function(x, lambda_vec, alphaX, betaX, gammaX) {
    #Params
    num_days = length(x)
    loglike = 0
    #prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
    
    for (t in 2:num_days) {
      #lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
      inner_sum_xt = 0
      term1_alpha = exp(-alphaX * lambda_vec[t])
      term2_alpha = alphaX * lambda_vec[t]
      
      for (yt in 0:x[t]) {
        #Sum for all values of yt
        
        #Log likelihood
        zt = x[t] - yt
        inner_sum_xt = inner_sum_xt +
          term1_alpha * (term2_alpha) ^ yt * #(1/factorial(yt))
          PROBABILITY_ZT(zt, lambda_vec[t], betaX, gammaX)
      }
      
      loglike = loglike + log(inner_sum_xt)
    }
    
    print(loglike)
    return(loglike)
  }

#2. PROBABILITY OF ZT
PROBABILITY_ZT <- function(zt, lambda_t, betaX, gammaX, max_nt = 5) {
  
  'Probability of Zt'
  
  #Initialise
  prob_zt = 0
  
  for (nt in 0:max_nt){
    #prob_zt = prob_zt + dpois(nt, betaX*lambda_t)*dpois(zt, gammaX*nt)
    prob_zt = prob_zt + poisson_density(nt, betaX*lambda_t)*poisson_density(zt, gammaX*nt)
  }
  
  return(prob_zt)
}

#LAMBDA FUNCTION
get_lambda <- function(epidemic_data, shape_gamma = 6, scale_gamma = 1){
  
  '#Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between '
  #Parameters
  num_days = length(epidemic_data)
  lambda_vec = vector("numeric", num_days)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Lambda -> days of the infection
  for (t in 1:num_days){
    lambda_vec[t] = sum(epidemic_data[1:(t-1)]*rev(prob_infect[1:(t-1)]))
  }
  
  return(lambda_vec)
}

#POISSON DENSITY
poisson_density <- function(x, rate){
  
  return((rate^x)*exp(-rate))
}
