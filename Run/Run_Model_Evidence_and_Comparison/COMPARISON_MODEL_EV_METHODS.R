#***************************
#COMPARE MODEL EVIDENCE; HARMONIC MEAN VS IMPORTANCE SAMPLING
#***************************

#PLOT
par(mfrow = c(2,1))

#1.COMPARE IS VS HARMONIC MEAN
#BASELINE MODEL
data_type = 'SSEB'; 
BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_base,
                                              harmonic_mean = list_hm_log_ev_base),
                      plot_title = paste0('Baseline Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))

#SSEB MODEL
BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_sseb,
                                              harmonic_mean = list_hm_log_ev_sseb),
                      plot_title = paste0('SSNB Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))


#SSNB MODEL
BOX_PLOT_MULT_RESULTS(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_ssnb,
                                              harmonic_mean = list_hm_log_ev_ssnb),
                      plot_title = paste0('SSNB Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))



#***************************
#1. IMPORTANCE SAMPLING COMPARISON
#***************************
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

