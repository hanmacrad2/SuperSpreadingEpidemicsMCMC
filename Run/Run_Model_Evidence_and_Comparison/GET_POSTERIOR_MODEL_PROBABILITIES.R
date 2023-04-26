
#***********************
# 5. POSTERIOR MODEL PROBABILITIES
#**********************
par(mfrow = c(1,1))

#RESULTS
data_type = 'SSEB'

#***********************
#MODEL 1: BASE
#***********************
vec_post_probs_hm_base = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_base,
                                                                    mod2 = list_is_log_ev_sseb, mod3 = list_is_log_ev_ssnb))

mean(vec_post_probs_hm_base); sd(vec_post_probs_hm_base)
#PLOT_BAYES_FACTORS(vec_post_probs_hm)

#***********************
#MODEL 2: SSEB
#***********************
vec_post_probs_hm_sseb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_sseb,
                                                                     mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_ssnb))
mean(vec_post_probs_hm_sseb); sd(vec_post_probs_hm_sseb)
#PLOT_BAYES_FACTORS(vec_post_probs_hm_sseb)

#***********************
#MODEL 3: SSNB
#***********************
vec_post_probs_ssnb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssnb,
                                                                     mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_sseb))
mean(vec_post_probs_ssnb); sd(vec_post_probs_ssnb)
#PLOT_BAYES_FACTORS(vec_post_probs_hm3)

#PLOT
BOX_PLOT_MULT_RESULTS(list_vec_results = list(SSEB = vec_post_probs_hm_sseb, BASE = vec_post_probs_hm_base,
                                              SSNB = vec_post_probs_ssnb),
                                              data_type = data_type)
