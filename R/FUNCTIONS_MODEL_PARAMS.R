#****************************************************************************
#* FUNCTIONS; MODELS & PARAMETERS

#********************************************
#1. MODELS
#*******************************************
GET_FLAGS_MODELS <-function(BASELINE = FALSE, SSE = FALSE, SSI = FALSE , SSEB = FALSE, SSIB = FALSE){
  
  if(BASELINE){
    FLAGS_MODELS = list(Baseline = TRUE, SSE = FALSE, SSI = FALSE,
                        SSEB = FALSE, SSIB = FALSE) 
  } else if (SSE){
    FLAGS_MODELS = list(Baseline = FALSE, SSE = TRUE, SSI = FALSE,
                        SSEB = FALSE, SSIB = FALSE) 
  } else if (SSI) {
    
    FLAGS_MODELS = list(Baseline = FALSE, SSE = FALSE, SSI = TRUE,
                        SSEB = FALSE, SSIB = FALSE) 
    
  } else if (SSEB) {
    
    FLAGS_MODELS = list(Baseline = FALSE, SSE = FALSE, SSI = FALSE,
                        SSEB = TRUE, SSIB = FALSE) 
    
  } else if (SSIB) {
    
    FLAGS_MODELS = list(Baseline = FALSE, SSE = FALSE, SSI = FALSE,
                        SSEB = FALSE, SSIB = TRUE) 
    
  }
  
  return(FLAGS_MODELS)
}

#***************************************************
#2. PARAMETERS
#***************************************************
GET_PARAM <- function(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE) {
  
  'Return parameter label'
  
  if(r0){
    FLAG_PARAM = list(r0 = TRUE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE)
    
  } else if (k){
    FLAG_PARAM = list(r0 = FALSE, k = TRUE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE)
    
  } else if (alpha) {
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = TRUE,
                      beta = FALSE, a = FALSE, b = FALSE, c = FALSE)
    
  } else if (beta){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = TRUE, a = FALSE, b = FALSE)
  } else if (a){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = TRUE, b = FALSE, c = FALSE)
  } else if (b){
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = TRUE, c = FALSE)
  } else if (c) {
    FLAG_PARAM = list(r0 = FALSE, k = FALSE, alpha = FALSE,
                      beta = FALSE, a = FALSE, b = FALSE, c = TRUE)
  }
  
  return(FLAG_PARAM)
}

#***************************************************
#3. PARAMETERS LABELS
#*******************************************************
GET_PARAM_LABEL <- function(FLAG_PARAM, model) { #model or FLAG_MODELS
  
  #PARAM
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  ylab_inf = 'Estimated posterior mean of '
  
  if(FLAG_PARAM$r0){
    list_labels = list(lab = expression(paste('R'[0])), 
                       xlab = expression(paste('True R'[0])), 
                       ylab = expression(paste('Estimated posterior mean of R'[0])), 
                       main_inf =  bquote(paste(italic(R[0]) ~ " - " ~ .(model), ' Model')), #. Estimated Posteriors & Prior')),
                       main_trace =  bquote(paste(italic(R[0]), " MCMC Trace")), #.(model)~ "model")),
                       main_hist = bquote(paste(italic(R[0]), " Posterior. ")),
                       main_hist_prior = bquote(paste(italic(R[0]), " Posterior. Prior: Exponential(1)")),
                       main_mean_sim = bquote(paste(italic(R[0]), " Cumulative mean. Simulated = 2.0")), #value 
                       main_mean0 = bquote(paste(italic(R[0]), " Cumulative mean - ", .(model)~ "model")),
                       legend_posterior = expression(paste('SSIB Model - Estimated Posteriors of R'[0], ' N = 1000')))
    
  } else if (FLAG_PARAM$alpha) {
    list_labels = list(lab = expression(paste(italic(alpha))), 
                       xlab = bquote(paste("True " ~ italic(alpha))),
                       ylab = bquote(paste(.(ylab_inf) ~ italic(alpha))),
                       main_inf = bquote(paste(italic(alpha) ~ " - " ~ .(model))),
                       main_trace =  bquote(paste(italic(alpha), " MCMC Trace")), #, .(model)~ "model")),
                       main_hist = bquote(paste(italic(alpha), " Posterior. ", .(model)~ "model")),
                       main_hist_prior = bquote(paste(italic(alpha), " Posterior. Prior: Beta(2,2)")),
                       main_mean_sim = bquote(paste(italic(alpha), " Cumulative mean. Simulated = 0.5")),
                       main_mean2 = bquote(paste(italic(alpha), " Cumulative mean - ", .(model)~ "model")),
                       legend_posterior = expression(paste("Estimated Posteriors of  ", alpha, " . N = 1000")))
    
  } else if (FLAG_PARAM$beta){
    list_labels = list(lab = expression(paste(italic(beta))), 
                       xlab = bquote(paste("True " ~ italic(beta))),
                       ylab = bquote(paste(.(ylab_inf) ~ italic(beta))),
                       main_inf =  bquote(paste(italic(beta) ~ " - " ~ .(model))),
                       main_trace =  bquote(paste(italic(beta), " MCMC Trace")), #.(model)~ "model")),
                       main_hist = bquote(paste(italic(beta), " Posterior - ", .(model)~ "model")),
                       main_hist_prior = bquote(paste(italic(beta), " Posterior. Prior: 1 + Gamma(3,3)")),
                       main_mean_sim = bquote(paste(italic(beta), " Cumulative mean. Simulated = 10")),
                       main_mean2 = bquote(paste(italic(beta), " Cumulative mean - ", .(model)~ "model")),
                       legend_posterior = expression(paste("Estimated Posteriors of  ", beta, " . N = 1000")))
  } else {
    
    list_labels = list(lab = bquote(paste(.(param))), 
                       xlab = paste0('True ', param),
                       ylab = paste0('Estimated posterior mean of ', param),
                       main_inf =  bquote(paste(.(param) ~ " - " ~ .(model))),
                       main_trace =  bquote(paste(.(param), " MCMC Trace")), #, .(model)~ "model")), #" Trace - ", .(model)~ "model")),
                       main_hist =  bquote(paste(.(param), " Posterior - ", .(model)~ "model")),
                       main_mean2 =  bquote(paste(.(param), " Cumulative mean - ", .(model)~ "model")),
                       legend_posterior = paste0('Estimated Posteriors of ', param, '. N = 1000'))
    
    list_labels = GET_ADDITIONAL_TITLES(FLAG_PARAM, list_labels)
  } 
  
  return(list_labels)
}

#***************************************************
#4. ADDITIONAL PARAMETERS LABELS
#***********************************************
GET_ADDITIONAL_TITLES <- function(FLAG_PARAM, list_labels){
  
  param = names(FLAG_PARAM)[which(unlist(FLAG_PARAM))]
  
  if(FLAG_PARAM$k){
    
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: Exponential(5)"))
    list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated = 0.1"))
    
  } else if (FLAG_PARAM$a) {
    
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: Beta(2,2)"))
    list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated = 0.5"))
    
  } else if (FLAG_PARAM$b) {
    
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: 1 + Gamma(3,3)"))
    list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated value = 10"))
    
  } else if (FLAG_PARAM$c){
    
    list_labels$main_hist_prior = bquote(paste(.(param), " Posterior. Prior: 1 + Gamma(3,3)"))
    list_labels$main_mean_sim = bquote(paste(.(param), " Cumulative mean. Simulated = 10"))
  }
  
  return(list_labels)
}

#***************************************************
#5. MODEL COLOURS
#***************************************************
GET_MODEL_COLORS <- function(){
  
  MODEL_COLORS <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red')
  
  return(MODEL_COLORS)
}

#***************************************************
#6. PRIORS
#***************************************************
GET_PRIOR <- function(EXP = FALSE, GAMMA = FALSE, UNIF = FALSE,
                       BETA = FALSE, GAMMA_B = FALSE) {
  
  if(EXP){
    FLAG_PRIOR = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE,
                      BETA = FALSE, GAMMA_B = FALSE)
    
  } else if (GAMMA){
    FLAG_PRIOR = list(EXP = FALSE, GAMMA = TRUE, UNIF = FALSE,
                      BETA = FALSE, GAMMA_B = FALSE)
    
  } else if (UNIF) {
    FLAG_PRIOR = list(EXP = FALSE, GAMMA = FALSE, UNIF = TRUE,
                      BETA = FALSE, GAMMA_B = FALSE)
    
  } else if (BETA){
    FLAG_PRIOR = list(EXP = FALSE, GAMMA = FALSE, UNIF = FALSE,
                      BETA = TRUE, GAMMA_B = TRUE)
  } else if (GAMMA_B){
    FLAG_PRIOR = list(EXP = FALSE, GAMMA = FALSE, UNIF = FALSE,
                      BETA = FALSE, GAMMA_B = FALSE)
  } 
  
  return(FLAG_PRIOR)
}