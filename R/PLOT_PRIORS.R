#*****************************
#* PLOT PRIORS - OVERLAPPING PRIOR ON FIGURES                       
#*****************************

PLOT_PRIOR_DIST <- function(FLAG_PARAM, xlimits = c(0,3), ylimits = c(0,15),
                            alpha = 0.2){ #0.4
  
  #PRIORS
  x_min = xlimits[1] #min(mcmc_vec) 
  x_max = xlimits[2] #max(mcmc_vec)
  x = seq(from = x_min, to = x_max, length = 5000)
              
  if(FLAG_PARAM$r0){
    x_min = 0; x_max = 10
    x = seq(from = x_min, to = x_max, length = 5000)
    y = dexp(x, 1)
    
  } else if (FLAG_PARAM$k){
    
    rate_k = GET_LIST_PRIORS_SSE()$k[1]
    print(paste0('rate_k', rate_k))
    
    y = dexp(x, rate = rate_k)
    
  } else if (FLAG_PARAM$alpha | FLAG_PARAM$a){
    
    x_min = 0; x_max = 1
    x = seq(from = x_min, to = x_max, length = 5000)
    y = dbeta(x, shape1 = 2, shape2 = 2)
    
  } else if (FLAG_PARAM$beta | FLAG_PARAM$b){ 
    
    x_min = 0 #min(mcmc_vec) 
    x_max = 20 # max(mcmc_vec)
    x = seq(from = x_min, to = x_max, length = 5000)
    y = dgamma(x-1, shape = 3, scale = 3)
  }
  
  #PLOT
  col_grey = rgb(0.5, 0.5, 0.5, alpha = alpha) 
  
  lines(x, y, type = 'l', lty = 1, col = col_grey)
  polygon(c(x, rev(x)), c(y, rep(0, length(y))),
          col = col_grey, border = NA)

}


#PRIOR TITLE
GET_PRIOR_TITLE <-function(FLAG_PARAM){
  
  #PRIORS
  if(FLAG_PARAM$r0){
    prior_title =  'Prior: Exponential(1),'
  } else if (FLAG_PARAM$k){
    prior_title =  'Prior: Exponential(5),'
  } else if (FLAG_PARAM$alpha | FLAG_PARAM$a){
    prior_title =  'Prior: Beta(2, 2),'
  } else if (FLAG_PARAM$beta | FLAG_PARAM$b){
    prior_title =  'Prior: 1 + Gamma(3, 3), '
  }
  
  return(prior_title)
}

#*********************************************************
#* PLOT PRIORS - THESIS                  
#*********************************************************
PRIOR_DIST_FIGURE <- function(PRIORS = list(EXP = FALSE, EXP_K = FALSE, EXP_K_5 = TRUE,
                                            GAMMA = FALSE, UNIF = FALSE, BETA_ALPHA = FALSE, BETA_A = TRUE,
                                            GAMMA_B = FALSE, GAMMA_BETA = FALSE),
                              PLOT_OVERLAP = TRUE, cex = 2.65){
  
  
  #PRIORS
  if(PRIORS$EXP){
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    y = dexp(x, 1)
    x_label = expression(R[0])
    prior_title =  bquote("Prior; " ~ p[R[0]] ~ "= Exponential(1)")
    #prior_title =  paste0('Exponential(1) Prior')
    
  } else if (PRIORS$EXP_K_5) {
    rate = 5
    x_min = 0; x_max = 1.5
    x = seq(from = x_min, to = x_max, length = 500)
    y = dexp(x, rate)
    x_label = expression(k)
    prior_title =  bquote("Prior; " ~ p[k] ~ "= Exponential(5)")
    
  } else if (PRIORS$EXP_K){
    rate = 1
    x_min = 0; x_max = 5
    x = seq(from = x_min, to = x_max, length = 500)
    y = dexp(x, rate)
    x_label = expression(k)
    prior_title =  bquote("Prior; " ~ p[k] ~ "= Exponential(1)")
    
  } else if (PRIORS$GAMMA){
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    prior_title =  paste0('Gamma(1, 5) Prior')
    y = dgamma(x, shape = 1, scale = 5)
    
  } else if (PRIORS$UNIF){
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    prior_title =  paste0('Uniform(0, 10) Prior')
    y = dunif(x, 0, 10)
    
  } else if (PRIORS$BETA_ALPHA){
    
    x_min = 0; x_max = 1
    x = seq(from = x_min, to = x_max, length = 500)
    #prior_title =  paste0('Beta(2, 2) Prior')
    #x_label = expression(beta)
    prior_title =  bquote("Prior; " ~ p[alpha] ~ "= Beta(2, 2)")
    y = dbeta(x, 2, 2)
    
  } else if (PRIORS$BETA_A) {
    
    x_min = 0; x_max = 1
    x = seq(from = x_min, to = x_max, length = 500)
    prior_title =  bquote("Prior; " ~ p[a] ~ "= Beta(2, 2)")
    y = dbeta(x, 2, 2)
    
  } else if (PRIORS$GAMMA_BETA){
    x_min = 0; x_max = 30
    x = seq(from = x_min, to = x_max, length = 500)
    x_label = expression(beta)
    prior_title =  bquote("Prior; " ~ p[beta] ~ "= 1 + Gamma(3, 3)")
    y = dgamma(x-1, shape = 3, scale = 3)
    #prior_title =  paste0('1 + Gamma(3, 3) Prior')
    
  } else if (PRIORS$GAMMA_B){
    
    x_min = 0; x_max = 25
    x = seq(from = x_min, to = x_max, length = 500)
    x_label = expression(b)
    prior_title =  bquote("Prior; " ~ p[b] ~ "= 1 + Gamma(3, 3)")
    y = dgamma(x-1, shape = 3, scale = 3)
  }
  
  #PLOT
  par(mfrow = c(2,3))
  par(mar=c(5.2, 5.2, 5.2, 5.2), xpd=TRUE)
  
  plot(x, y, type = 'l', lwd = 2, #col = 'orange',
       main = prior_title, 
       #xlab = x_label, 
       ylab = 'Density',
       cex.lab=cex, cex.axis=cex-0.3, cex.main= cex, cex.sub=cex-0.3) 
  
}

#TITLE
GET_PRIOR_TITLE_FIGURE <-function(PRIORS){
  
  #PRIORS
  if(PRIORS$EXP){
    prior_title =  'Exponential(1) Prior used'
  } else if (PRIORS$GAMMA){
    prior_title =  'Gamma(1, 5) Prior used'
  } else if (PRIORS$UNIF){
    prior_title =  'Uniform(0, 10) Prior used'
  } else if (PRIORS$BETA){
    prior_title =  'Beta(2, 2) Prior used'
  } else if (PRIORS$GAMMA_B){
    prior_title =  '1 + Gamma(3, 3) Prior used' 
  }
}

#***********
#PRIOR CIs
library(stats)

GET_PRIOR_CI <- function(FLAG_PARAM){
  
  if(FLAG_PARAM$r0){
    lower_bound <- 0.025; 
    upper_bound <- 3.688
    prior_mean = 1
    x0 <- lower_bound - 0.1
   
  } else if (FLAG_PARAM$k) {
    lower_bound <- 0.051; upper_bound <- 0.693
    
  } else if (FLAG_PARAM$alpha){
    lower_bound <- 0.035; upper_bound <- 0.965
    prior_mean = 0.5
    x0 <- x0 - 0.02
  } else if (FLAG_PARAM$a) {
    lower_bound <- 0.035; upper_bound <- 0.965
    prior_mean = 0.5
    
  } else if (FLAG_PARAM$beta) {
    prior_ci = GET_GAMMA_CI()
    #x0 <- x_lim[1] - 0.25
  }
  else if (FLAG_PARAM$b){
    prior_ci = GET_GAMMA_CI()
    #x0 <- x_lim[1] - 0.5
  } 
  
  prior_ci = c(lower_bound, upper_bound, prior_mean)
  return(prior_ci)
}

GET_GAMMA_CI <- function(shape = 3, scale = 3){
  
  # Calculate the lower and upper bounds of the 95% confidence interval
  lower_bound <- 1 + qgamma(0.025, shape, scale=scale)
  upper_bound <- 1 + qgamma(0.975, shape, scale=scale)
  mean_gamma = 1 + shape*scale
  cat("mean_gamma= ", mean_gamma, ", 95% Confidence Interval:", lower_bound, "->", upper_bound, "\n")
  gamma_ci = c(lower_bound, upper_bound, mean_gamma)
  
  return(gamma_ci)
}


