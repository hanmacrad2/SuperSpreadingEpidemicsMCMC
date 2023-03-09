#***************************************
#1. GET Q
#***************************************
GET_LOG_Q_PROPOSAL_UNI_VAR <- function(mcmc_samples, epidemic_data, 
                                       n_samples, num_dims = 1) {}

#MUTLI DIM PROPOSAL
GET_LOG_Q_PROPOSAL_MULTI_DIM <- function(mcmc_samples, epidemic_data,  
                                         n_samples) {}

#************************************************************
#* 2. GET P_HATS (log)
#************************************************************
GET_LOG_P_HAT <- function(mcmc_samples, epidemic_data, n_samples = 10000) {}

#*********************
#* 3. GET BAYES FACTORS
#*********************
GET_LOG_BAYES_FACTORS <- function (list_log_phat_mod1, list_log_phat_mod2){}
  
#*****************************************
#* 4.  GET POSTERIOR MODEL PROBABILITIES
#******************************************
GET_POSTERIOR_MODEL_PROBS <- function(num_models = 3, 
                                     probs_models = list(prob1 = 0.25, probh2 = 0.5, prob3 = 0.25), 
                                     log_phats = list(mod1 = mod1,
                                                      mod2 = mod2, mod3 = mod3)){}

#*************************
# 5. LOAD MCMC + GET POSTERIOR MODEL PROBABILITIES
#*****************************************************

LOAD_MCMC_GET_P_HAT <- function(epidemic_data, OUTPUT_FOLDER, run = 1, n_repeats = 100,
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                    SSIB = FALSE, SSIC = FALSE))