#**********************************************
#
# MODEL EVIDENCE ESTIMATES VIA IMPORTANCE SAMPLING
#
#***********************
run = 1 
n_repeats = 5

EPI_DATA = data_sseb
data_type = 'SSEB Data' #,  k= 0.1, P(k) ~ exp(1), P(k) ~ exp(0.1)' #BASE data' #SSE data'  

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
#2. SSE
#*************************
model_ev_sse = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                        SSIB = FALSE, SSIR = FALSE))

mean(model_ev_sse)
sd(model_ev_sse)

#*************************
#3. SSI
#*************************
model_ev_ssi = LOAD_MCMC_GET_SSIR_MODEL_EV(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats)

mean(model_ev_ssi)
sd(model_ev_ssi)

#EXP_01
run = 'exp_01'
model_ev_ssi_exp01 = LOAD_MCMC_GET_SSIR_MODEL_EV(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats)

mean(model_ev_ssi_exp01)
sd(model_ev_ssi)

#*************************
#4. SSE-B
#*************************
model_ev_sseb = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSE = FALSE,
                                                                        SSIB = FALSE, SSI = FALSE))
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
#PLOT ALL MODEL EVIDENCE
#*************************
par(mfrow = c(2,1))
BOX_PLOT_MODEL_EV(list_vec_results = list(BASE = model_ev_base, SSE = model_ev_sse, 
                                          SSI = model_ev_ssi, SSEB = model_ev_sseb,
                                         SSIB = model_ev_ssib), data_type = data_type)

BOX_PLOT_MODEL_EV(list_vec_results = list(BASE = model_ev_base, SSE = model_ev_sse, 
                                          SSI = model_ev_ssi), data_type = data_type)


BOX_PLOT_MODEL_EV(list_vec_results = list(SSI_exp1 = model_ev_ssi, SSI_exp01 = model_ev_ssi_exp01),
                  data_type = data_type)
