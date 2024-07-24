#*******************************************************
#*
# 1. PRIORS - DEFINE LISTS
#*
#*****************************************************************
#BASELINE
GET_LIST_PRIORS_BASELINE <- function() {
  
  LIST_PRIORS_BASELINE = list(exp = c(1,0), #exp
                         gamma = c(1, 5),
                         r0_unif = c(0,10),
                         unif =  c(0,10)) #exp(1), exp(0.1)
  
  return(LIST_PRIORS_BASELINE)
}

#SSE
GET_LIST_PRIORS_SSE <- function() {
  
  LIST_PRIORS_SSE = list(r0 = c(1,0), #exp
                         r0_gamma = c(1, 5),
                         r0_unif = c(0,10),
                    k =  c(5,0),
                    k_unif = c(0, 2)) #c(1,0) #exp(1), exp(0.1)
  
  return(LIST_PRIORS_SSE)
}

#SSI
GET_LIST_PRIORS_SSI <- function(){
  
  LIST_PRIORS_SSI = list(r0 = c(1,0), #r0 = list(exp = c(1,0), gamma = c(1,5))
                         gamma = c(1,5),
                         r0_unif = c(0,10),
                    k =  c(5,0),
                    k_unif = c(0, 2)) #c(1,0) exp
  
  return(LIST_PRIORS_SSI)
}

#SSE-B
GET_LIST_PRIORS_SSEB <- function() {
  
  LIST_PRIORS_SSEB = list(r0 = c(1,0), #r0; exp dist 
                          r0_unif = c(0,10),
                     alpha =  c(2,2), #alpha; beta dist alpha =  c(1,2)
                     alpha_unif = c(0,1),
                     beta = c(3,3),
                     beta_unif = c(1,40)) #beta; gamma dist [2,1], [8,1]
  
  return(LIST_PRIORS_SSEB)
}

#SSI-B
GET_LIST_PRIORS_SSIB <- function(){
  
  LIST_PRIORS_SSIB = list(r0 = c(1,0),    #exp dist 
                          r0_unif = c(0,10),
                     a =  c(2,2), #beta dist c(2,2)
                     a_unif = c(0,1),
                     b = c(3,3),
                     b_unif = c(1,40)) #gamma dist (2,1)
  
  return(LIST_PRIORS_SSIB)
}

#GET PRIORS USED
GET_PRIORS_USED <- function(){
  
  PRIORS_USED = 
    list(BASELINE = 
           list(r0 = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE)),
         SSE = 
           list(r0 = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE),
                k =  list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE)),
         SSI = 
           list(r0 = list(EXP = TRUE, GAMMA = FALSE, UNIF = FALSE),
                k =  list(EXP = TRUE, UNIF = FALSE)),
         SSEB =
           list(r0 = list(EXP = TRUE, UNIF = FALSE),
                alpha = list(BETA = TRUE, EXP = FALSE, UNIF = FALSE),
                beta = list(GAMMA = TRUE, EXP = FALSE, UNIF = FALSE),
                gamma = FALSE),
         SSIB = 
           list(r0 = list(EXP = TRUE, UNIF = FALSE),
                a = list(BETA = TRUE, EXP = FALSE, GAMMA = FALSE, UNIF = FALSE),
                b = list(EXP = FALSE, GAMMA = TRUE, UNIF = FALSE)), 
         c = FALSE)
  
  return(PRIORS_USED)
}

#SET UNIFORM PRIORS
SET_UNIFORM_PRIORS <- function(PRIORS_USED){
  
  'Set Uniform Priors'
  
  #1. BASELINE MODEL
  PRIORS_USED$BASELINE$r0$UNIF = TRUE
  PRIORS_USED$BASELINE$r0$EXP = FALSE
  
  #2. SSE MODEL
  PRIORS_USED$SSE$r0$UNIF = TRUE
  PRIORS_USED$SSE$r0$EXP = FALSE
  
  PRIORS_USED$SSE$k$UNIF = TRUE
  PRIORS_USED$SSE$k$EXP = FALSE
  
  #3. SSI
  PRIORS_USED$SSI$r0$UNIF = TRUE
  PRIORS_USED$SSI$r0$EXP = FALSE
  
  PRIORS_USED$SSI$k$UNIF = TRUE
  PRIORS_USED$SSI$k$EXP = FALSE
  
  #4. SSEB MODEL
  PRIORS_USED$SSEB$r0$UNIF = TRUE
  PRIORS_USED$SSEB$r0$EXP = FALSE
  
  PRIORS_USED$SSEB$alpha$UNIF = TRUE
  PRIORS_USED$SSEB$alpha$BETA = FALSE
  
  PRIORS_USED$SSEB$beta$UNIF = TRUE
  PRIORS_USED$SSEB$beta$GAMMA = FALSE
  
  #5. SSIB MODEL
  PRIORS_USED$SSIB$r0$UNIF = TRUE
  PRIORS_USED$SSIB$r0$EXP = FALSE
  
  PRIORS_USED$SSIB$a$UNIF = TRUE
  PRIORS_USED$SSIB$a$BETA = FALSE
  
  PRIORS_USED$SSIB$b$UNIF = TRUE
  PRIORS_USED$SSIB$b$GAMMA = FALSE
  
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
