
#***********************
# 5. POSTERIOR MODEL PROBABILITIES
#**********************
par(mfrow = c(1,1))

#RESULTS
data_type = 'SSEB'
data_type = 'BASELINE'
model_ev_method = 'IS'

#*************************
#1. IMPORTANCE SAMPLING DATA
#*************************

#***********************
#MODEL 1: BASE
#***********************
vec_post_probs_is_base1 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_base1,
                                                                    mod2 = list_is_log_ev_sseb1, mod3 = list_is_log_ev_ssnb1))

mean(vec_post_probs_is_base1); sd(vec_post_probs_is_base1)
#PLOT_BAYES_FACTORS(vec_post_probs_hm)

#***********************
#MODEL 2: SSEB
#***********************
vec_post_probs_is_sseb1 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_sseb1,
                                                                     mod2 = list_is_log_ev_base1, mod3 = list_is_log_ev_ssnb1))
mean(vec_post_probs_is_sseb1); sd(vec_post_probs_is_sseb1)
#PLOT_BAYES_FACTORS(vec_post_probs_hm_sseb)

#***********************
#MODEL 3: SSNB
#***********************
vec_post_probs_is_ssnb1 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssnb,
                                                                     mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_sseb))
mean(vec_post_probs_is_ssnb1); sd(vec_post_probs_is_ssnb1)
#PLOT_BAYES_FACTORS(vec_post_probs_hm3)

#PLOT POSTERIOR PROBS
model_ev_method = 'IS'
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(BASE = vec_post_probs_is_base1, SSEB = vec_post_probs_is_sseb1,
                                              SSNB = vec_post_probs_is_ssnb1),
                                              data_type = data_type, model_ev_method = model_ev_method)

#Data
plot.ts(data_baseline, ylab = 'Infection count', main = paste0(data_type, ' Simulated data'))

#*************************
#2. HARMONIC MEAN DATA
#*************************
model_ev_method = 'HM'

#***********************
#MODEL 1: BASE
#***********************
vec_post_probs_hm_base = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_hm_log_ev_base,
                                                                         mod2 = list_hm_log_ev_sseb, mod3 = list_hm_log_ev_ssnb))

mean(vec_post_probs_hm_base); sd(vec_post_probs_hm_base)
#PLOT_BAYES_FACTORS(vec_post_probs_hm)

#***********************
#MODEL 2: SSEB
#***********************
vec_post_probs_hm_sseb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_hm_log_ev_sseb,
                                                                         mod2 = list_hm_log_ev_base, mod3 = list_hm_log_ev_ssnb))
mean(vec_post_probs_hm_sseb); sd(vec_post_probs_hm_sseb)
#PLOT_BAYES_FACTORS(vec_post_probs_hm_sseb)

#***********************
#MODEL 3: SSNB
#***********************
vec_post_probs_hm_ssnb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_hm_log_ev_ssnb,
                                                                      mod2 = list_hm_log_ev_base, mod3 = list_hm_log_ev_sseb))
mean(vec_post_probs_hm_ssnb); sd(vec_post_probs_hm_ssnb)
#PLOT_BAYES_FACTORS(vec_post_probs_hm3)

#PLOT POSTERIOR PROBS
model_ev_method = 'HM'
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSEB = vec_post_probs_hm_sseb, BASE = vec_post_probs_hm_base,
                                              SSNB = vec_post_probs_hm_ssnb),
                      data_type = data_type, model_ev_method = model_ev_method)
#Data
plot.ts(data_sseb, ylab = 'Infection count', main = paste0(data_type, ' Simulated data'))
