% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSE_MCMC.R
\name{SSE_MCMC_ADAPTIVE}
\alias{SSE_MCMC_ADAPTIVE}
\title{MCMC adaptive algorithm for Super-Spreading Events epidemic model}
\usage{
SSE_MCMC_ADAPTIVE(
  epidemic_data,
  mcmc_inputs = list(n_mcmc = 5e+05, mod_start_points = list(m1 = 0.8, m2 = 0.05, m3 =
    10), alpha_star = 0.4, thinning_factor = 10),
  priors_list = list(alpha_prior_exp = c(1, 0), beta_prior_ga = c(10, 2/100),
    beta_prior_exp = c(0.1, 0), gamma_prior_ga = c(10, 1), gamma_prior_exp = c(0.1, 0)),
  FLAGS_LIST = list(ADAPTIVE = TRUE, ABG_TRANSFORM = TRUE, PRIOR = TRUE, BETA_PRIOR_GA
    = FALSE, GAMMA_PRIOR_GA = FALSE, THIN = TRUE)
)
}
\arguments{
\item{epidemic_data}{data from the epidemic, namely daily infection counts}

\item{mcmc_inputs}{A list of mcmc specifications including
\itemize{
  \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
  \item \code{"mod_start_points"} - Model parameter starting points; where the mcmc algorithm begins sampling from
  \item \code{"alpha_star"} - target acceptance rate; used in the adaptive algorithm
  \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
}}

\item{priors_list}{A list of prior parameters used
\itemize{
  \item \code{"alpha_prior_exp"} - rate of exponential prior on alpha, default rate of 1
  \item \code{"beta_prior_exp"} - rate of exponential prior on beta, default rate of 1
  \item \code{"beta_prior_ga"} -  shape and scale of gamma distribution prior on beta, defaults Ga(10, 0.02)
  \item \code{"gamma_prior_exp"}  - rate of exponential prior on gamma, default rate of 0.1
  \item \code{"gamma_prior_ga"} -  shape and scale of gamma distribution prior on c, defaults Ga(10, 1)
}}

\item{FLAGS_LIST}{A list of Boolean variables for switching on/off certain functionality
\itemize{
  \item \code{"ADAPTIVE"} - Adaptive MCMC algorithm implemented if TRUE
  \item \code{"ABG_TRANSFORM"} - Transform update step proposed as part of MCMC algorithm if TRUE **
  \item \code{"PRIOR"}  - Apply prior distributions to model parameters
  \item \code{"BETA_PRIOR_GA"}  - A Gamma prior on beta if TRUE, otherwise exponential
  \item \code{"GAMMA_PRIOR_GA"}  - A Gamma prior on gamma if TRUE, otherwise exponential
  \item \code{"THIN"}  - Return a thinned MCMC sample if TRUE, reduced by a factor of \code{"thinning_factor"}
}}
}
\value{
MCMC samples, acceptance rates etc.
\itemize{
  \item \code{"alpha_vec"} - A vector containing MCMC samples of the SSI model parameter alpha. Vector size \code{"mcmc_vec_size = n_mcmc/thinning_factor"}
  \item \code{"beta_vec"} - A vector containing MCMC samples of the SSI model parameter beta. Vector size \code{"mcmc_vec_size"}
  \item \code{"gamma_vec"} - A vector containing MCMC samples of the SSI model parameter gamma. Vector size \code{"mcmc_vec_size"}
  \item \code{"r0_vec"} - A vector containing samples of the SSI model parameter r0 obtained via; r0 = alpha + beta*gamma. Vector size \code{"mcmc_vec_size"}
  \item \code{"log_like_vec"}  - A vector containing the log likelihood at each of the MCMC iterations, of size of \code{"mcmc_vec_size"}
  \item \code{"sigma"} - A list of vectors containing the adapted sigma for each of the model parameters at each iteration of the MCMC, if \code{"FLAGS_LIST$ADAPTIVE"} is TRUE. Otherwise sigma is a list of constant sigma values for each of the model parameters
  \item \code{"list_accept_rates"}  - A list of the MCMC acceptance rates (percentage) for each of the model parameters
  }
}
\description{
MCMC algorithm with Adaptation for obtaining samples from the parameters of a
super-spreading events (SSE) epidemic model
}
\examples{

mcmc_inputs = list(n_mcmc = 500000, mod_start_points = list(m1 = 0.8, m2 = 0.05, m3 = 10),
alpha_star = 0.4, thinning_factor = 10)

#START MCMC
mcmc_sse_output = SSE_MCMC_ADAPTIVE(epidemic_data, mcmc_inputs)

}
\author{
Hannah Craddock, Xavier Didelot, Simon Spencer
}
