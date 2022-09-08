#********************
#UTIL FUNCTIONS

#'Log-sum-exp of a vector
#'
#' Log-sum-exp trick applied to a vector to overcome the issue of under- or overflow.
#'  Often arise when normalizing vectors of log probabilities
#' 
#' @param vectorX A vector of quantities 
#' @return normalized exponentiated vector
#' @export LOG_SUM_EXP
#'
#' @author Hannah Craddock
#'
#' @examples
#'
#' vectorX = LOG_SUM_EXP(vectorX)
#'

#LOG EXP SUM
LOG_SUM_EXP <- function(vectorX){
  
  max_val = max(vectorX)
  
  out = max_val + log(sum(exp(vectorX - max_val)))
  
  return(out)
}


ll = c(-1000, -1000, -1000)
lse = LOG_SUM_EXP(ll)
ll_norm = exp(ll - lse) 
ll_norm