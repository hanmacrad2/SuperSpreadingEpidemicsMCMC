
#***********************
# 5. POSTERIOR MODEL PROBABILITIES
#**********************
par(mfrow = c(1,1))

#***********************
#MODEL 1: BASE
#***********************
vec_post_probs_hm = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_log_ev_base,
                                                                    mod2 = list_log_ev_sseb, mod3 = list_log_ev_ssnb))

mean(vec_post_probs_hm); sd(vec_post_probs_hm)
PLOT_BAYES_FACTORS(vec_post_probs_hm)

#***********************
#MODEL 2: SSEB
#***********************
vec_post_probs_hm2 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_log_ev_sseb,
                                                                     mod2 = list_log_ev_base, mod3 = list_log_ev_ssnb))
mean(vec_post_probs_hm2); sd(vec_post_probs_hm2)
PLOT_BAYES_FACTORS(vec_post_probs_hm2)

#***********************
#MODEL 3: SSNB
#***********************
vec_post_probs_hm3 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_log_ev_ssnb,
                                                                     mod2 = list_log_ev_base, mod3 = list_log_ev_sseb))
mean(vec_post_probs_hm3); sd(vec_post_probs_hm3)
PLOT_BAYES_FACTORS(vec_post_probs_hm3)
