#*****************************
#* PLOT PRIORS                       
#*****************************
PLOT_PRIOR_DIST <- function(PRIORS = list(EXP = FALSE, EXP_K = FALSE,
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
    
  } else if (PRIORS$EXP_K){
    x_min = 0; x_max = 5
    x = seq(from = x_min, to = x_max, length = 500)
    y = dexp(x, 1)
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
  
  if(PLOT_OVERLAP){
    
    par(mfrow = c(2,3))
    par(mar=c(5.2, 5.2, 5.2, 5.2), xpd=TRUE)
    
  } else {
    plot(x, y, type = 'l', lwd = 2, #col = 'orange',
         main = prior_title, 
         #xlab = x_label, 
         ylab = 'Density',
         cex.lab=cex, cex.axis=cex-0.3, cex.main= cex, cex.sub=cex-0.3) 
  }
  
}


#*****************************
#* PLOT PRIORS                       
#*****************************
PRIOR_DIST_FIGURE <- function(PRIORS = list(EXP = FALSE, EXP_K = FALSE,
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
    
  } else if (PRIORS$EXP_K){
    x_min = 0; x_max = 5
    x = seq(from = x_min, to = x_max, length = 500)
    y = dexp(x, 1)
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
  
  if(PLOT_OVERLAP){
    
    par(mfrow = c(2,3))
    par(mar=c(5.2, 5.2, 5.2, 5.2), xpd=TRUE)
    
  } else {
    plot(x, y, type = 'l', lwd = 2, #col = 'orange',
         main = prior_title, 
         #xlab = x_label, 
         ylab = 'Density',
         cex.lab=cex, cex.axis=cex-0.3, cex.main= cex, cex.sub=cex-0.3) 
  }

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
  } else if (PRIORS$beta | FLAG_PARAM$b){
    prior_title =  'Prior: 1 + Gamma(3, 3)'
  }
}

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