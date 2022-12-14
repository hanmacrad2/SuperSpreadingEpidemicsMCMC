% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSI_MCMC.R
\name{SSI_MCMC_ADAPTIVE}
\alias{SSI_MCMC_ADAPTIVE}
\title{MCMC adaptive algorithm for Super-Spreading Individuals epidemic model}
\usage{
SSI_MCMC_ADAPTIVE(
  data,
  mcmc_inputs = list(n_mcmc = 5e+05, mod_start_points = list(m1 = 0.72, m2 = 0.0038, m3
    = 22), alpha_star = 0.4, thinning_factor = 10),
  priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp =
    c(0.1, 0), c_prior_ga = c(10, 1), c_prior_exp = c(0.1, 0)),
  FLAGS_LIST = list(ADAPTIVE = TRUE, DATA_AUG = TRUE, BCA_TRANSFORM = TRUE, PRIOR =
    TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE, THIN = TRUE)
)
}
\arguments{
\item{data}{data from the epidemic, namely daily infection counts}

\item{mcmc_inputs}{A list of mcmc specifications including
\itemize{
  \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler
  \item \code{"mod_start_points"} - Model parameter starting points; where the mcmc algorithm begins sampling from
  \item \code{"alpha_star"} - target acceptance rate; used in the adaptive algorithm
  \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
}}

\item{priors_list}{A list of prior parameters used
\itemize{
  \item \code{"a_prior_exp"} - rate of exponential prior on a, default rate of 1
  \item \code{"b_prior_exp"} - rate of exponential prior on b, default rate of 1
  \item \code{"b_prior_ga"} -  shape and scale of gamma distribution prior on b, defaults Ga(10, 0.02)
  \item \code{"c_prior_exp"}  - rate of exponential prior on c, default rate of 0.1
  \item \code{"c_prior_ga"} -  shape and scale of gamma distribution prior on c, defaults Ga(10, 1)
}}

\item{FLAGS_LIST}{A list of Boolean variables for switching on/off certain functionality
\itemize{
  \item \code{"ADAPTIVE"} - Adaptive MCMC algorithm implemented if TRUE
  \item \code{"DATA_AUG"} - Data Augmentation implemented as part of the SSI model if TRUE
  \item \code{"BCA_TRANSFORM"} - Transform update step proposed as part of MCMC algorithm if TRUE **
  \item \code{"PRIOR"}  - Apply prior distributions to model parameters
  \item \code{"B_PRIOR_GAMMA"}  - A Gamma prior on b if TRUE, otherwise exponential
  \item \code{"C_PRIOR_GAMMA"}  - A Gamma prior on c if TRUE, otherwise exponential
  \item \code{"THIN"}  - Return a thinned MCMC sample if TRUE, reduced by a factor of \code{"thinning_factor"}
}}
}
\value{
MCMC samples, acceptance rates etc.
\itemize{
  \item \code{"a_vec"} - A vector containing MCMC samples of the SSI model parameter a. Vector size  mcmc_vec_size = \code{"n_mcmc/thinning_factor"}
  \item \code{"b_vec"} - A vector containing MCMC samples of the SSI model parameter b. Vector size \code{"mcmc_vec_size"}
  \item \code{"c_vec"} - A vector containing MCMC samples of the SSI model parameter c. Vector size \code{"mcmc_vec_size"}
  \item \code{"r0_vec"} - A vector containing samples of the SSI model parameter r0 obained via; r0 = a + b*c. Vector size \code{"mcmc_vec_size"}
  \item \code{"log_like_vec"}  - A vector containing the log likelihood at each of the mcmc iterations, of size of \code{"mcmc_vec_size"}
  \item \code{"sigma"} - A list of vectors containing the adapted sigma for each of the model parameters at each iteration of the MCMC, if \code{"FLAGS_LIST$ADAPTIVE"} is TRUE. Otherwise sigma is a list of constant sigma values for each of the model parameters
  \item \code{"list_accept_rates"}  - A list of the MCMC acceptance rates (%) for each of the model parameters
  \item \code{"non_ss"}  - The augmented epidemic data corresponding to the non super-spreaders (non-ss) for each iteration of the mcmc
  \item \code{"ss"}  - The augmented epidemic data corresponding to the super-spreaders (ss) for each iteration of the mcmc
  }
}
\description{
MCMC algorithm with Adaptation for obtaining samples from the parameters of a
super-spreading individuals (SSI) epidemic model. Includes; Adaptation, Data Augmentation, B-C & A-C transform
}
\examples{

mcmc_inputs = list(n_mcmc = 500000, mod_start_points = list(m1 = 0.72, m2 = 0.0038, m3 = 22),
alpha_star = 0.4, thinning_factor = 10)

#START MCMC
mcmc_ssi_output = SSI_MCMC_ADAPTIVE(epidemic_data, mcmc_inputs)

}
\author{
Hannah Craddock, Xavier Didelot, Simon Spencer
}
