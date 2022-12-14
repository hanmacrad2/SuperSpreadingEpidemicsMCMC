% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODEL_CRITICISM.R
\name{GET_P_VALUES_TOTAL}
\alias{GET_P_VALUES_TOTAL}
\title{Get and save the posterior predictive p value (ppp-value) for each summary statistic}
\usage{
GET_P_VALUES_TOTAL(root_folder, n_reps)
}
\arguments{
\item{root_folder}{root folder location in which to store the MCMC results for \code{"n_reps"}}

\item{n_reps}{Total number of repetitions of the MCMC sampler ran. Ppp-values calculated on the aggregate of all reps.}

\item{epidemic_data}{data from the epidemic, namely daily infection counts}
}
\value{
A Dataframe of ppp-values for all the summary statistics \code{"df_p_values"}. Also saves the dataframe in the \code{"root_foler"} location
}
\description{
Get and save the posterior predictive p value (ppp-value) for each summary statistic across all reps. The ppp-value of each summary statistic is calculated and stored together in a dataframe.
}
\examples{

GET_P_VALUES_TOTAL(root_folder, n_reps)

}
\author{
Hannah Craddock, Xavier Didelot, Simon Spencer
}
