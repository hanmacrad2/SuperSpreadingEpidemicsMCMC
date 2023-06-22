#RUN POSTERIOR PROBS II
post_probs_base = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base2,
                                                                                                                model_ev_ssnb2, model_ev_sseb2,
                                                                                                                model_ev_ssir2, model_ev_ssib2))

post_probs_ssnb = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_ssnb2,
                                                                                          model_ev_base2, model_ev_sseb2,
                                                                                          model_ev_ssir2, model_ev_ssib2))



post_probs_sseb = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_sseb2,
                                                                                          model_ev_base2, model_ev_ssnb2,
                                                                                          model_ev_ssir2, model_ev_ssib2))

post_probs_ssir = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_ssir2,
                                                                                          model_ev_base2, model_ev_ssnb2, 
                                                                                          model_ev_sseb2, model_ev_ssib2))

post_probs_ssib = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_ssib2,
                                                                                          model_ev_base2, model_ev_ssnb2, 
                                                                                          model_ev_sseb2, model_ev_ssir2))
#PLOT
#model_ev_method = 'IS'
data_type = 'SSI-B Data' #MOCK DATA' #'# (20% Missing)'
#data_type = 'NZ Waitemata 08/21 Subset I (3 models)'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results =
                           list(BASE = post_probs_base, SSI = post_probs_ssir,
                                SSNB = post_probs_ssnb,
                                SSEB = post_probs_sseb, SSIB = post_probs_ssib),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')
