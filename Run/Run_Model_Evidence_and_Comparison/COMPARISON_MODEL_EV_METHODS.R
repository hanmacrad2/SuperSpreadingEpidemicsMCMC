#COMPARE MODEL EVIDENCE; HARMONIC MEAN VS IMPORTANCE SAMPLING

#PLOT
par(mfrow = c(1,1))

BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_base,
                                              harmonic_mean = list_log_ev_base),
                      plot_title = 'Baseline Model Evidence. IS alg (2018) vs HM. 500 runs, N mcmc = 500k')

BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_sseb3,
                                              harmonic_mean = list_log_ev_sseb),
                      plot_title = 'SSEB Model Evidence. Importance sampling alg (2018) vs Harmonic Mean')


                                     