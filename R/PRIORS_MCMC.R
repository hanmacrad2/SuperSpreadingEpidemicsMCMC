#*******************************************************
#*
# 1. PRIORS - DEFINE LISTS
#*
#*****************************************************************
#BASELINE
GET_LIST_PRIORS_BASELINE <- function() {
  
  LIST_PRIORS_BASELINE = list(exp = c(1,0), #exp
                         gamma = c(1, 5),
                         unif =  c(0,10)) #exp(1), exp(0.1)
  
  return(LIST_PRIORS_BASELINE)
}

#SSE
GET_LIST_PRIORS_SSE <- function() {
  
  LIST_PRIORS_SSE = list(r0 = c(1,0), #exp
                         r0_gamma = c(1, 5),
                         r0_unif = c(0,10),
                    k =  c(1,0)) #exp(1), exp(0.1)
  
  return(LIST_PRIORS_SSE)
}

#SSI
GET_LIST_PRIORS_SSI <- function(){
  
  LIST_PRIORS_SSI = list(r0 = c(1,0), #r0 = list(exp = c(1,0), gamma = c(1,5))
                         gamma = c(1,5),
                    k =  c(1,0)) #exp
  
  return(LIST_PRIORS_SSI)
}

#SSE-B
GET_LIST_PRIORS_SSEB <- function() {
  
  LIST_PRIORS_SSEB = list(r0 = c(1,0), #r0; exp dist 
                     alpha =  c(2,2), #alpha; beta dist alpha =  c(1,2)
                     gamma = c(3,3)) #beta; gamma dist [2,1], [8,1]
  
  return(LIST_PRIORS_SSEB)
}

#SSI-B
GET_LIST_PRIORS_SSIB <- function(){
  
  LIST_PRIORS_SSIB = list(r0 = c(1,0),    #exp dist 
                     alpha =  c(2,2), #beta dist c(2,2)
                     b = c(3,3)) #gamma dist (2,1)
  
  return(LIST_PRIORS_SSIB)
}

#GET PRIORS USED
GET_PRIORS_USED <- function(){
  
  PRIORS_USED = 
    list(BASELINE = 
           list(r0 = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE)),
         SSE = 
           list(r0 = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE),
                k =  list(EXP = TRUE, GAMMA = FALSE)),
         SSI = 
           list(r0 = list(EXP = TRUE),
                k =  list(EXP = TRUE)),
         SSEB =
           list(r0 = list(EXP = TRUE),
                alpha = list(EXP = FALSE, BETA = TRUE),
                beta = list(EXP = FALSE, GAMMA = TRUE),
                gamma = FALSE),
         SSIB = 
           list(r0 = list(EXP = TRUE),
                alpha = list(EXP = FALSE, BETA = TRUE, GAMMA = FALSE),
                b = list(EXP = FALSE, GAMMA = TRUE)), 
         c = FALSE)
  
  return(PRIORS_USED)
}

#SET_PRIORS
SET_PRIORS <- function(){
  
  return(list(list_priors = list(
    priors_sse = GET_LIST_PRIORS_SSE(),
    priors_ssi = GET_LIST_PRIORS_SSI(),
    priors_sseb = GET_LIST_PRIORS_SSEB(),
    priors_ssib = GET_LIST_PRIORS_SSIB()
  ), PRIORS_USED = GET_PRIORS_USED()))
}
