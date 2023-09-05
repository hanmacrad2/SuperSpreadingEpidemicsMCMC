
# Super-spreading in Epidemics: A Modelling & Model Comparison Framework

This package contains the code to implement a modelling framework encompassing five models and inter-model comparison that aims to account for a range of different spreading dynamics. Each of the models are discrete-time, stochastic, branching-process models that aim to model the spreading of an epidemic through a susceptible population. The five-model framework encompasses a baseline model with regular-spreading, a model with individual reproduction number (SSI), a bimodal super-spreading events (SSEB) model, a bimodal super-spreading individuals (SSIB) model and a negative binomial model with over-dispersion to account for super-spreading (SSE). The basic reproductive rate (R0) is a key measure in estimating the ability of a new pathogen to spread, and is a parameter in each of the models. A Bayesian approach is used to infer model parameters by means of Markov chain Monte Carlo (MCMC). The framework is adaptive and generalizable to a wide range of diseases and is applied to both simulated data and real data (Covid-19, SARS outbreak (2003)). Model comparison is carried out between the models to determine which model best fits a particular dataset.

The package also contains other util functions such as plotting functions to plot the mcmc samples and results


## Authors

- [@hanmacrad2](https://www.github.com/hanmacrad2)
- [@drsimonspencer](https://github.com/drsimonspencer)
- [@xavierdidelot](https://github.com/xavierdidelot)

