#************************************
#GET POSTERIOR MODEL PROBABILITIES
#*************************************

#****************************************************
# DATASET 1
#***************************************************

#************************
# PART 1: FIVE MODELS :D
#************************
post_probs_base = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base,
                                                                                                                  model_ev_ssnb, model_ev_sseb,
                                                                                                                model_ev_ssir, model_ev_ssib))

post_probs_ssnb = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_ssnb,
                                                                                            model_ev_base, model_ev_sseb,
                                                                                          model_ev_ssir, model_ev_ssib))



post_probs_sseb = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_sseb,
                                                                                            model_ev_base, model_ev_ssnb,
                                                                                          model_ev_ssir, model_ev_ssib))

post_probs_ssir = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_ssir,
                                                                                          model_ev_base, model_ev_ssnb, 
                                                                                          model_ev_sseb, model_ev_ssib))

post_probs_ssib = GET_AGG_POSTERIOR_PROBABILITES(num_models = 5, list_log_mod_evid = list(model_ev_ssib,
                                                                                          model_ev_base, model_ev_ssnb, 
                                                                                          model_ev_sseb, model_ev_ssir))
#PLOT
#model_ev_method = 'IS'
data_type = 'SSE Data'
#data_type = 'NZ Waitemata 08/21 Subset I (3 models)'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results =
                           list(BASE = post_probs_base, SSE = post_probs_ssnb, SSI = post_probs_ssir,
                                                 SSEB = post_probs_sseb, SSIB = post_probs_ssib),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')


#************************
# PART 1: TWO MODELS
#************************
post_probs_base2 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 2, FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base,
                                                                                                                 model_ev_ssir))

post_probs_ssir2 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 2, list_log_mod_evid = list(model_ev_ssir,
                                                                                           model_ev_base))

#PLOT
model_ev_method = 'IS'
data_type = 'NZ Waitemata 08/21 Subset I (2 models)'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSI = post_probs_ssir2, BASE = post_probs_base2),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')

#************************
# PART 2: THREE MODELS
#************************
post_probs_base3 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base,
                                                                                                model_ev_ssnb,
                                                                                                model_ev_sseb))

post_probs_ssnb3 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, list_log_mod_evid = list(model_ev_ssnb,
                                                                          model_ev_ssir,
                                                                          model_ev_ssib))



post_probs_ssib3 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, list_log_mod_evid = list(model_ev_ssib,
                                                                                           model_ev_ssnb,
                                                                                           model_ev_ssir))

post_probs_ssir3 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, list_log_mod_evid = list(model_ev_ssir,
                                                                          model_ev_ssnb,
                                                                          model_ev_ssib))

#PLOT
#model_ev_method = 'IS'
#data_type = 'Simulated Baseline'# (20% Missing)'
#data_type = 'NZ Waitemata 08/21 Subset I (3 models)'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSIB = post_probs_ssib3,
                                                SSI = post_probs_ssir3, SSNB = post_probs_ssnb3),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')



#************************
# PART 3: FOUR MODELS
#************************

post_probs_base = GET_AGG_POSTERIOR_PROBABILITES(FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base,
                                                                                 model_ev_sseb,
                                                                                 model_ev_ssnb,
                                                                                 model_ev_ssir))

post_probs_sseb = GET_AGG_POSTERIOR_PROBABILITES(list_log_mod_evid = list(model_ev_sseb,
                                                                          model_ev_base,
                                                                          model_ev_ssnb,
                                                                          model_ev_ssir))

post_probs_ssnb = GET_AGG_POSTERIOR_PROBABILITES(list_log_mod_evid = list(model_ev_ssnb,
                                                                          model_ev_base,
                                                                          model_ev_sseb,
                                                                          model_ev_ssir))

post_probs_ssir = GET_AGG_POSTERIOR_PROBABILITES(list_log_mod_evid = list(model_ev_ssir,
                                                                          model_ev_base,
                                                                          model_ev_sseb,
                                                                          model_ev_ssnb))

#PLOT
model_ev_method = 'IS'
data_type = 'NZ Waitemata 08/21 Subset I'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSEB = post_probs_sseb, BASE = post_probs_base,
                                                 SSNB = post_probs_ssnb, SSI = post_probs_ssir),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')


#****************************************************
# DATASET 2
#***************************************************

#************************
# PART 1: TWO MODELS
#************************
post_probs_base22 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 2, FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base2,
                                                                                                                 model_ev_ssir2))

post_probs_ssir22 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 2, list_log_mod_evid = list(model_ev_ssir2,
                                                                                           model_ev_base2))

#PLOT
model_ev_method = 'IS'
data_type = 'NZ Waitemata 08/21 Subset I (2 models)'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSI = post_probs_ssir2, BASE = post_probs_base2),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')

#************************
# PART 2: THREE MODELS
#************************
post_probs_base32 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, FLAG_BASELINE = TRUE, list_log_mod_evid = list(model_ev_base2,
                                                                                                                 model_ev_ssnb2,
                                                                                                                 model_ev_ssir2))

post_probs_ssnb32 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, list_log_mod_evid = list(model_ev_ssnb2,
                                                                                           model_ev_base2,
                                                                                           model_ev_ssir2))

post_probs_ssir32 = GET_AGG_POSTERIOR_PROBABILITES(num_models = 3, list_log_mod_evid = list(model_ev_ssir2,
                                                                                           model_ev_base2,
                                                                                           model_ev_ssnb2))

#PLOT
model_ev_method = 'IS'
data_type = 'NZ Waitemata 08/21 Subset II (3 models)'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSNB = post_probs_ssnb32,
                                                 SSI = post_probs_ssir32, BASE = post_probs_base32),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')

#************************
# PART 4: FOUR MODELS
#************************

post_probs_base = GET_AGG_POSTERIOR_PROBABILITES(4, FLAG_BASELINE = TRUE,
                                                  list_log_mod_evid = list(model_ev_base,
                                                                           model_ev_sseb,
                                                                           model_ev_ssnb,
                                                                           model_ev_ssib))

post_probs_sseb = GET_AGG_POSTERIOR_PROBABILITES(4, list_log_mod_evid = list(model_ev_sseb,
                                                                          model_ev_base,
                                                                          model_ev_ssnb,
                                                                          model_ev_ssib))

post_probs_ssnb = GET_AGG_POSTERIOR_PROBABILITES(4, list_log_mod_evid = list(model_ev_ssnb,
                                                                          model_ev_base,
                                                                          model_ev_sseb,
                                                                          model_ev_ssir))

post_probs_ssib = GET_AGG_POSTERIOR_PROBABILITES(4, list_log_mod_evid = list(model_ev_ssir,
                                                                          model_ev_base,
                                                                          model_ev_sseb,
                                                                          model_ev_ssnb))

#PLOT
model_ev_method = 'IS'
data_type = 'Mock Data' #NZ Waitemata 08/21 Subset II'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSEB = post_probs_sseb, BASE = post_probs_base,
                                                 SSNB = post_probs_ssnb, SSIB = post_probs_ssib),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')

