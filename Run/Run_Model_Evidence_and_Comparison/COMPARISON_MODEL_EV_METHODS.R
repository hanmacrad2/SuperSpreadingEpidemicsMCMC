#COMPARE MODEL EVIDENCE; HARMONIC MEAN VS IMPORTANCE SAMPLING

#PLOT
par(mfrow = c(1,1))

#PARAMETERS(RENAME)

#1. COMPARE IS VS HARMONIC MEAN
#BASELINE MODEL
BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_base,
                                              harmonic_mean = list_log_ev_base),
                      plot_title = 'Baseline Model Evidence. IS alg (2018) vs HM. 500 runs, N mcmc = 500k')

#SSNB MODEL
BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_ssnb,
                                              harmonic_mean = list_log_ev_ssnb),
                      plot_title = 'SSNB Model Evidence. IS alg (2018) vs HM. 500 runs, N mcmc = 500k')

#SSEB MODEL
BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_sseb,
                                              harmonic_mean = list_log_ev_sseb),
                      plot_title = 'SSEB Model Evidence. Importance sampling alg (2018) vs Harmonic Mean')


#2. COMPARE SAMPLES RUNS

#***************************
#1. IMPORTANCE SAMPLING COMPARISON
BOX_PLOT_MULT_RESULTS(list_vec_results = list(model_ev_500k = ests_phat_base,
                                              model_ev_50k = ests_phat_base_100),
                      plot_title = 'Important Sampling (2018)
                      N mcmc = 500k, 500 runs  vs  N mcmc = 50k, 100 runs.  Baseline model')

BOX_PLOT_MULT_RESULTS(list_vec_results = list(model_ev_500k = ests_phat_ssnb,
                                              model_ev_50k = ests_phat_ssnb_100),
                      plot_title = 'Important Sampling (2018)
                     N mcmc = 500k, 500 runs  vs  N mcmc = 50k, 100 runs.  SSNB model')

#***************************
#2. HARMONIC MEAN COMPARISON
BOX_PLOT_MULT_RESULTS(list_vec_results = list(model_ev_500k = list_log_ev_base,
                                              model_ev_50k = list_log_ev_base_100),
                      plot_title = 'Harmonic Mean
                      N mcmc = 500k, 500 runs  vs  N mcmc = 50k, 100 runs. Baseline model')

BOX_PLOT_MULT_RESULTS(list_vec_results = list(model_ev_500k = list_log_ev_ssnb,
                                              model_ev_50k = list_log_ev_ssnb_100),
                      plot_title = 'Harmonic Mean
                     N mcmc = 500k, 500 runs  vs  N mcmc = 50k, 100 runs. SSNB model')

