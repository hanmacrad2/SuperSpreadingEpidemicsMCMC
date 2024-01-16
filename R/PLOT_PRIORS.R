#*****************************
#* PLOT PRIORS                       
#*****************************

PLOT_PRIOR_DIST <- function(FLAG_PARAM, mcmc_vec, limits){
  
  #PRIORS
  x_min = limits[1] #min(mcmc_vec) 
  x_max = limits[2] #max(mcmc_vec)
  x = seq(from = x_min, to = x_max, length = 5000)
              
  if(FLAG_PARAM$r0){
    #x_min = 0.9; x_max = 3.0
    #x = seq(from = x_min, to = x_max, length = 5000)
    y = dexp(x, 1)
    
  } else if (FLAG_PARAM$k){
    
    y = dexp(x, 1)
    
  } else if (FLAG_PARAM$alpha | FLAG_PARAM$a){
    
    #x_min = 0; x_max = 1
    #x = seq(from = x_min, to = x_max, length = 5000)
    y = dbeta(x, shape1 = 2, shape2 = 2)
    
  } else if (FLAG_PARAM$beta | FLAG_PARAM$b){
    
    #x_min = min(mcmc_vec) 
    #x_max = max(mcmc_vec)
    #x = seq(from = x_min, to = x_max, length = 5000)
    y = dgamma(x-1, shape = 3, scale = 3)
  }
  
  #PLOT
  col_grey = rgb(0.5, 0.5, 0.5, alpha = 0.4)
  lines(x, y, type = 'l', lty = 1, col = col_grey,  ylim = c(0, max(y)))
  
  #Filled area
  polygon(c(x, rev(x)), c(y, rep(min(y), length(y))),
          col = col_grey, border = NA)
}


#PRIOR TITLE
GET_PRIOR_TITLE <-function(FLAG_PARAM){
  
  #PRIORS
  if(FLAG_PARAM$r0){
    prior_title =  'Prior: Exponential(1)'
  } else if (FLAG_PARAM$k){
    prior_title =  'Prior: Exponential(1)'
  } else if (FLAG_PARAM$alpha | FLAG_PARAM$a){
    prior_title =  'Prior: Beta(2, 2)'
  } else if (FLAG_PARAM$beta | FLAG_PARAM$b){
    prior_title =  'Prior: 1 + Gamma(3, 3)'
  }
  
  return(prior_title)
}
