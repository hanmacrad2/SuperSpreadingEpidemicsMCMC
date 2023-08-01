#**********************************************
#
# MODEL EVIDENCE ESTIMATES VIA IMPORTANCE SAMPLING
#
#***********************
run = 1 
n_repeats = 5

EPI_DATA = data_baseline

#*************************
#1. BASELINE
#*************************
model_ev_base = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = FALSE))
mean(model_ev_base)
sd(model_ev_base) 
#PLOT_MODEL_EV_RESULTS(model_ev_base, model_type = '. Baseline w/ exp(1) prior on R0')

#*************************
#2. SSNB
#*************************
model_ev_ssnb = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                        SSIB = FALSE, SSIR = FALSE))

mean(model_ev_ssnb)
sd(model_ev_ssnb)

#*************************
#3. SSEB
#*************************
model_ev_sseb = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = FALSE))
#PLOT_MODEL_EV_RESULTS(model_ev_sseb)
mean(model_ev_sseb)
sd(model_ev_sseb)

#*************************
#5. SSIB
#*************************
model_ev_ssib = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                             FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                 SSIB = TRUE, SSIR = FALSE))
mean(model_ev_ssib)
sd(model_ev_ssib)

#*************************
#4. SSIR
#*************************
model_ev_ssir = LOAD_MCMC_GET_SSIR_MODEL_EV(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats)

mean(model_ev_ssir)
sd(model_ev_ssir)

#*************************
#PLOT ALL MODEL EVIDENCE
#*************************
par(mfrow = c(2,1))
data_type = 'Baseline'  
BOX_PLOT_MODEL_EV(list_vec_results = list(BASE = model_ev_base, SSE = model_ev_ssnb, 
                                          SSI = model_ev_ssir, SSEB = model_ev_sseb,
                                         SSIB = model_ev_ssib), data_type = data_type)


#***************************************************************************
#*
#* GAMMA PRIORS
#* 
#**************************************************************************

#BASELINE + GAMMA PRIOR 
run = '1_ga'
model_ev_base_gp = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                                                    SSIB = FALSE, SSIR = FALSE))
mean(model_ev_base_gp)# , na.rm = TRUE)
sd(model_ev_base_gp) #, na.rm = TRUE)
PLOT_MODEL_EV_RESULTS(model_ev_base_gp, model_type = '. Baseline w/ gamma(2, mean = R0) prior on R0')

#*************************
#SSIB + GAMMA PRIOR 
#*************************
run = '2_ga_prior'
model_ev_ssib2 = LOAD_MCMC_GET_MODEL_EVIDENCE(data_ssib2, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                             FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                 SSIB = TRUE, SSIR = FALSE))
mean(model_ev_ssib2)
sd(model_ev_ssib2)
PLOT_MODEL_EV_RESULTS(model_ev_ssib2, model_type = 'SSIB; prior(a) ~ gamma(mean ~ a)',
                      data_type = 'SSIB')

#*************************
#SSNB + GAMMA PRIOR ON K
#*************************
run = 'gamma_prior_k'; 
model_ev_ssnb_gak = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub2, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                       FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                           SSIB = FALSE, SSIR = FALSE))

PLOT_MODEL_EV_RESULTS(model_ev_ssnb_gak)

#PLOT
data_type = 'Waitemata days 10-20'
par(mfrow = c(1,1))
BOX_PLOT_MODEL_EV(list_vec_results = list(ssnb_k_exp_prior = model_ev_ssnb,
                                          ssnb_k_ga_prior = model_ev_ssnb_gak),
                              title = 'SSNB Model k priors compared exp(1) vs Ga(0.001, rt = 0.001). ',
                              data_type = data_type, model = '') 

#MODEL EVIDENCE
mcmc_samples = cbind(mcmc_output_ssir$ssir_params_matrix, mcmc_output_ssir$eta_matrix)
list_mod_ev = GET_LOG_MODEL_EVIDENCE(mcmc_samples, EPI_DATA, FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                                 SSIB = FALSE, SSIR = TRUE))

#SSIR

#MATRIX
mat_ssir = matrix(0, nrow = 5, ncol = 10)

BOX_PLOT_RESULTS(mat_ssir3, title = 'Model Evidence ',
                 model_type = 'SSIR ',
                 data_type = 'Mock Data, 6 days')


#*****************************
#*
#LOAD MODEL EVIDENCE ESTIMATES
#*
#*****************************
LOAD_MODEL_EVIDENCE <- function(model_type, run, OUTER_FOLDER){
  

  print(model_type)
  CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  list_model_ev = readRDS(file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
  print(list_model_ev)
  return(list_model_ev)
}

#LOAD 
run = 2
model_ev_base2 = LOAD_MODEL_EVIDENCE('baseline', run, OUTER_FOLDER)
model_ev_sseb2 = LOAD_MODEL_EVIDENCE('sseb', run, OUTER_FOLDER)
model_ev_ssnb2 = LOAD_MODEL_EVIDENCE('ssnb', run, OUTER_FOLDER)
model_ev_ssir2 = LOAD_MODEL_EVIDENCE('ssir', run, OUTER_FOLDER)

#SD OF RESULTS
sd(model_ev_ssir2, na.rm = TRUE) 
