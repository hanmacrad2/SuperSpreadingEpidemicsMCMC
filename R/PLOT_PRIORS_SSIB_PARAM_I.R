#PLOT PRIORS - SSIB PARAM I

#*****************************
#* PLOT PRIORS                       
#*****************************

PLOT_PRIOR_DIST_SSIB_I <- function(FLAG_PARAM, mcmc_vec, limits){
  
  #PRIORS
  x_min = limits[1] #min(mcmc_vec) 
  x_max = limits[2] #max(mcmc_vec)
  x = seq(from = x_min, to = x_max, length = 5000)
  
  if(FLAG_PARAM$a){
    x_min = 0.9; x_max = 3.0
    x = seq(from = x_min, to = x_max, length = 5000)
    y = dexp(x, 1)
    
  } else if (FLAG_PARAM$b){
    
    y = dexp(x, 1)
    
  } else if (FLAG_PARAM$c){ 
    
    x_min = 0 #min(mcmc_vec) 
    x_max = 20 # max(mcmc_vec)
    x = seq(from = x_min, to = x_max, length = 5000)
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
GET_PRIOR_TITLE_SSIB_I <-function(FLAG_PARAM){
  
  #PRIORS
  if(FLAG_PARAM$a | FLAG_PARAM$b){
    prior_title =  'Prior: Exponential(1)'
  } else if (FLAG_PARAM$c){
    prior_title =  'Prior: 1 + Gamma(3, 3)'
  }
  
  return(prior_title)
}

#MCMC HISTOGRAM + #PLOT PRIOR TITLE
PLOT_MCMC_HIST_SSIB_I <- function (mcmc_vec, FLAGS_MODELS, FLAG_PARAM, MODEL_COLOR,
                            cex = 1.6, xlim = c(1.0, 3.0)){
  
  model = names(FLAGS_MODELS)[which(unlist(FLAGS_MODELS))]
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  list_labels = GET_PARAM_LABEL(FLAG_PARAM, model)
  prior = GET_PRIOR_TITLE_SSIB_I(FLAG_PARAM)
  
  hist(mcmc_vec, freq = FALSE, breaks = 200,
       xlab = list_labels$lab,
       xlim = xlim,
       border = MODEL_COLOR,
       col = MODEL_COLOR, 
       main = list_labels$main_hist_prior,
       cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
}