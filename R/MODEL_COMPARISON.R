#MODEL COMPARISON

#*********************
#* MODEL COMPARISON BY RATIO OF MODEL EVIDENCES
#***********
MODEL_EVIDENCE <- function(loglike_vec){
  
  'Model evidence via log-sum-exp trick'
  
  m = max(loglike_vec, na.rm = TRUE)
  log_model_ev = m + log(mean(exp(loglike_vec - m)))
  
  return(log_model_ev)
}

get_bayes_factor <- function(loglike_vec1, loglike_vec2){
  
  'Get Bayes factor via ratio of the model evidence'
  
  bayes_factor = MODEL_EVIDENCE(loglike_vec1)/MODEL_EVIDENCE(loglike_vec2)
  
  return(bayes_factor)
}

#'Harmonic mean of Log likelihood
#'
#' Harmonic mean of Log likelihood. Used in conjunction with \code{"get_bayes_factor_harmonic_means"} to get the Bayes Factor between two competing models
#'
#' @param log_like_vec Log Likelihood at each iteration of the MCMC algorithm
#'
#' @return \code{"harmonic_mean"}  harmonic mean of log likelihood vector
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' harmonic_mean = get_harmonic_mean(log_like_vec)
#'
#'
#Get harmonic mean
get_harmonic_mean <- function(likelihood_vec){

  'Returns harmonic mean of the vector'

  return(1/(mean(1/likelihood_vec)))

}

#'Bayes Factor via Harmonic mean
#'
#'Bayes Factor of two competing models obtained via the ratio of their two harmonic means
#'
#' @param harmonic_mean1 harmonic mean obtained from model 1
#' @param harmonic_mean2 harmonic mean obtained from model 2
#'
#' @return \code{"bayes_factor"} bayes factor of two models up for comparison. Obtained via their harmonic means
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' bayes_factor = get_bayes_factor_hm(harmonic_mean1, harmonic_mean2)
#'
#'
#Get harmonic mean
get_bayes_factor_hm <- function(harmonic_mean1, harmonic_mean2){

  'Returns bayes factor; ratio of two harmonic means'
  return(harmonic_mean1/harmonic_mean2)

}

#'Bayes Factor via Harmonic mean
#'
#'Bayes Factor of two competing models obtained via the ratio of their two harmonic means
#'
#' @param log_like_vec1 Log Likelihood of model 1 at each iteration of the MCMC algorithm
#' @param log_like_vec2 Log Likelihood of model 2 at each iteration of the MCMC algorithm
#'
#' @return \code{"bayes_factor"} bayes factor of two models up for comparison. Obtained via their harmonic means
#'
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' bayes_factor = get_bayes_factor_hm(loglike_vec1, loglike_vec2)
#'
#'
#Get harmonic mean
get_bayes_factor_harmonic_means <- function(loglike_vec1, loglike_vec2){

  'Returns bayes factor; ratio of two harmonic means'
  harmonic_mean1 = get_harmonic_mean(loglike_vec1)
  harmonic_mean2 = get_harmonic_mean(loglike_vec2)

  bayes_factor = harmonic_mean1/harmonic_mean2

  return(bayes_factor)

}

# model_comparison_hm <- function(){
#
#
# }
