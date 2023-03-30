#COMPARE MODEL EVIDENCE; HARMONIC MEAN VS IMPORTANCE SAMPLING

BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_base2,
                                              harmonic_mean = list_log_ev_base),
                      plot_title = 'Baseline Model Evidence. Importance sampling alg (2018) vs Harmonic Mean')

BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = ests_phat_sseb,
                                              harmonic_mean = list_log_ev_sseb),
                      plot_title = 'SSEB Model Evidence. Importance sampling alg (2018) vs Harmonic Mean')


                                     