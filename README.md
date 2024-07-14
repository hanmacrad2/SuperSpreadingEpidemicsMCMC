
# Super-spreading in Epidemics: A Bayesian Modelling Framework with Multi-Model Comparison

This package contains the code to implement a bayesian modelling framework of epidemic transmission that encompasses five distinct models. The package also contains the code to implement multi-model comparison between the candidate models using a method of importance sampling to estimate the model evidence of each model. 

## Installation

You can install the package from CRAN using:

```r
install.packages("SuperSpreadingEpidemicsMCMC")
```
Or you can install the development version from GitHub:

```
install.packages("devtools")
devtools::install_github("hanmacrad2/SuperSpreadingEpidemicsMCMC")
```

## Background

The transmission dynamics of an epidemic are rarely homogeneous. Super-spreading events (SSEs) and super-spreading individuals (SSIs) are two examples of heterogeneous transmissibility. Inference of super-spreading is commonly carried out on secondary case data, the distribution of which is known as the offspring distribution $Z$. However, this data is often unavailable. The negative binomial distribution is often chosen to quantify super-spreading as it allows for over-dispersion in $Z$, but this model fails to capture extreme super-spreading as it is unimodal. Addressing this, we introduce a five-model framework fit to incidence time-series, data that is more commonly recorded and available. The framework consists of five discrete-time, stochastic, branching-process models of epidemic spread through a susceptible population. It incorporates novel, bimodal super-spreading models, is disease-agnostic and implemented as an R package. The five models include a Baseline (Poisson) model of homogeneous transmission, two super-spreading events models; a negative binomial model (SSE) and a novel bimodal model (SSEB) and two super-spreading individuals models, again a negative binomial model (SSI) and a novel bimodal model (SSIB). A Bayesian framework is used to infer model parameters using Markov Chain Monte Carlo (MCMC). Multi-model comparison, or model selection, is conducted by comparing the models' posterior probabilities, in which an importance sampling estimator is used to estimate the models' marginal likelihoods. The estimator is selected for its consistency and lower variance compared to the alternative harmonic mean estimator. Model selection, when applied to simulated data from each model, identifies the correct model for the majority of simulations and accurately infers the true parameters, including the basic reproduction number $R_0$. We apply our methods of model inference, comparison and criticism to incidence data from the Covid-19 pandemic caused by SARS-CoV-2 and the 2003 SARS outbreak. Model selection consistently identifies the same model and mechanism for a given disease, even when using different time series. Our estimates of the dispersion parameter $k$ of the negative binomial models align with those from studies that use secondary case data. Quantifying the contribution of super-spreading to disease transmission has important implications for infectious disease management and control. The novel modelling framework developed has potential to be a valuable tool for public health. 

This package contains the code to implement a modelling framework encompassing five models and multi-model comparison that aims to account for a range of different spreading dynamics. Each of the models are discrete-time, stochastic, branching-process models that aim to model the spreading of an epidemic through a susceptible population. The five-model framework encompasses a baseline model with regular-spreading, a model with individual reproduction number (SSI), a bimodal super-spreading events (SSEB) model, a bimodal super-spreading individuals (SSIB) model and a negative binomial model with over-dispersion to account for super-spreading (SSE). The basic reproductive rate (R0) is a key measure in estimating the ability of a new pathogen to spread, and is a parameter in each of the models. A Bayesian approach is used to infer model parameters by means of Markov chain Monte Carlo (MCMC). The framework is adaptive and generalizable to a wide range of diseases and is applied to both simulated data and real data (Covid-19, SARS outbreak (2003)). Model comparison is carried out between the models to determine which model best fits a particular dataset.

The package also contains other util functions such as plotting functions to plot the mcmc samples and results

## Example
library(SuperSpreadingEpidemicsMCMC)
data_baseline <- SIMULATE_EPI_BASELINE(r0 = 2.0)
data_sse <- SIMULATE_EPI_SSE(r0 = 2.0, k = 0.2)
data_ssi <- SIMULATE_EPI_SSI(r0 = 2.0, k = 0.2)
data_sseb <- SIMULATE_EPI_SSE(r0 = 2.0, alpha = 0.5, beta = 10)
data_ssib <- SIMULATE_EPI_SSIB(r0 = 2.0, a = 0.5, b = 10)
mcmc_baseline <- MCMC_INFER_BASELINE(data_baseline)
mcmc_sse <- MCMC_INFER_SSE(data_sse)
mcmc_ssi <- MCMC_INFER_SSI(data_ssi)
mcmc_sseb <- MCMC_INFER_SSEB(data_sseb)
mcmc_ssib <- MCMC_INFER_SSIB(data_ssib)


## Authors

- [@hanmacrad2](https://www.github.com/hanmacrad2)
- [@drsimonspencer](https://github.com/drsimonspencer)
- [@xavierdidelot](https://github.com/xavierdidelot)

