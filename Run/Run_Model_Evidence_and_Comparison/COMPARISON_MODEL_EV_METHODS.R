#***************************
#COMPARE MODEL EVIDENCE; HARMONIC MEAN VS IMPORTANCE SAMPLING
#***************************

#PLOT
par(mfrow = c(3,1))

#1.COMPARE IS VS HARMONIC MEAN
#BASELINE MODEL
data_type = 'SSEB'; 
data_type = 'NZ Waitemata Subset'
BOX_PLOT_MODEL_EV(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_base,
                                              harmonic_mean = list_hm_log_ev_base),
                  data_type = data_type, model = 'BASE')
                      #plot_title = paste0('Baseline Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))

#SSEB MODEL
BOX_PLOT_MODEL_EV(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_sseb,
                                              harmonic_mean = list_hm_log_ev_sseb),
                  data_type = data_type, model = 'SSEB')
                      #plot_title = paste0('SSEB Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))


#SSNB MODEL
BOX_PLOT_MODEL_EV(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_ssnb,
                                              harmonic_mean = list_hm_log_ev_ssnb),
                  data_type = data_type, model = 'SSNB')
                      #plot_title = paste0('SSNB Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))

#SSIR MODEL
BOX_PLOT_MODEL_EV(list_vec_results = list(importance_sampling_2018 = list_is_log_ev_ssir,
                                          harmonic_mean = list_hm_log_ev_ssir),
                  data_type = data_type, model = 'SSIR')
#plot_title = paste0('SSNB Model Evidence. ', data_type, ' data. IS alg (2018) vs HM. ', n_repeats, ' runs, N mcmc = 50k'))


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

