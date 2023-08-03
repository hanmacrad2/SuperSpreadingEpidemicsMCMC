#**********************************************
#
# MODEL EVIDENCE ESTIMATES VIA IMPORTANCE SAMPLING
#
#***********************
run = 2 
n_repeats = 5

#EPI_DATA = data_base
data_type = 'SSI'  

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
#5. SSI-B
#*************************
model_ev_ssib = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                             FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                 SSIB = TRUE, SSIR = FALSE))
mean(model_ev_ssib)
sd(model_ev_ssib)

#*************************
#4. SSI
#*************************
model_ev_ssir = LOAD_MCMC_GET_SSIR_MODEL_EV(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats)

mean(model_ev_ssir)
sd(model_ev_ssir)

#*************************
#PLOT ALL MODEL EVIDENCE
#*************************
par(mfrow = c(2,1))
BOX_PLOT_MODEL_EV(list_vec_results = list(BASE = model_ev_base, SSE = model_ev_ssnb, 
                                          SSI = model_ev_ssir, SSEB = model_ev_sseb,
                                         SSIB = model_ev_ssib), data_type = data_type)
