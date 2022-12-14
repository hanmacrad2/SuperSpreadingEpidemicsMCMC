% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODEL_CRITICISM.R
\name{GET_P_VALUE}
\alias{GET_P_VALUE}
\title{P value derived by comparing the true summary statistics of the true epidemic data with that of the summary statistics of the simulated data}
\usage{
GET_P_VALUE(column_true_val, column_summary_stat)
}
\arguments{
\item{column_true_val}{Column of the true summary stat values of the epidemic data}

\item{column_summary_stat}{Column of summary stats of the simulated data}
}
\value{
p-value
}
\description{
Returns the p value
}
\examples{

p_value = GET_P_VALUE(column_true_val, column_summary_stat)


}
\author{
Hannah Craddock, Xavier Didelot, Simon Spencer
}
