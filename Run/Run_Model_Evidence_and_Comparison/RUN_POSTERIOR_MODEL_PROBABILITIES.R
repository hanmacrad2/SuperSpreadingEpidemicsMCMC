library(SuperSpreadingEpidemicsMCMC)
#***********************
# 5. POSTERIOR MODEL PROBABILITIES
#**********************
par(mfrow = c(2,1))
par(mfrow = c(1,1))

#RESULTS
data_type = 'SSEB x30 different simulated'
data_type = 'BASELINE'
data_type = 'Waitemata NZ, August 2021 (T = 10)'
model_ev_method = 'IS'

#*************************
#1. IMPORTANCE SAMPLING DATA
#*************************

#***********************
#MODEL 1: BASE
#***********************
vec_post_probs_is_base = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_base,
                                                                    mod2 = list_is_log_ev_sseb, mod3 = list_is_log_ev_ssnb))

mean(vec_post_probs_is_base); sd(vec_post_probs_is_base)
#PLOT_BAYES_FACTORS(vec_post_probs_hm)

#***********************
#MODEL 2: SSEB
#***********************
vec_post_probs_is_sseb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_sseb,
                                                                     mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_ssnb),
                                                 probs_models = list(prob1 = 0.25, prob2 = 0.5, prob3 = 0.25))

mean(vec_post_probs_is_sseb); sd(vec_post_probs_is_sseb)
#PLOT_BAYES_FACTORS(vec_post_probs_hm_sseb)

#***********************
#MODEL 3: SSNB
#***********************
vec_post_probs_is_ssnb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssnb,
                                                                     mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_sseb),
                                                probs_models = list(prob1 = 0.25, prob2 = 0.5, prob3 = 0.25))
mean(vec_post_probs_is_ssnb); sd(vec_post_probs_is_ssnb)
#PLOT_BAYES_FACTORS(vec_post_probs_hm3)

#PLOT POSTERIOR PROBS
model_ev_method = 'IS'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(BASE = vec_post_probs_is_base, SSE_B = vec_post_probs_is_sseb,
                                              SSE_NB = vec_post_probs_is_ssnb),
                                              data_type = data_type, model_ev_method = '') #IS Model Evidence.')

#Data
plot.ts(data_wait_08_21_sub1, ylab = 'Infection count', main = paste0(data_type)) #, ', real data'))

#********************
#COMPARISON: TWO MODELS
#********************
post_probs_base_2 = GET_AGG_POSTERIOR_PROB_II(list_log_mod_evid = list(mod1 = list_is_log_ev_base,
                                                                       mod2 = list_is_log_ev_ssnb))

post_probs_ssnb_2 = GET_AGG_POSTERIOR_PROB_II(list_log_mod_evid = list(mod1 = list_is_log_ev_ssnb,
                                                   mod2 = list_is_log_ev_base))

BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(BASE = post_probs_base_2,
                                                 SSNB = post_probs_ssnb_2),
                         data_type = data_type, model_ev_method = '', titleX = 'Comparison of x2 models; ')

#Data
plot.ts(data_wait_08_21, ylab = 'Infection count', main = paste0(data_type, ', real data'))

#BF
#bf = GET_BAYES_FACTORS(list_is_log_ev_ssnb, list_is_log_ev_base)
#boxplot(bf)

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


#************************
# FOUR MODELS
#************************
#BASE
vec_post_probs_is_base = GET_AGG_POSTERIOR_PROB(FLAG_BASELINE = TRUE, list_log_mod_evid = list(mod1 = list_is_log_ev_base,
                                                                                               mod2 = list_is_log_ev_sseb, mod3 = list_is_log_ev_ssnb,
                                                                                               mod4 = list_is_log_ev_ssir))

#SSEB
vec_post_probs_is_sseb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_sseb,
                                                                         mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_ssnb,
                                                                         mod4 = list_is_log_ev_ssir))

#SSNB
vec_post_probs_is_ssnb = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssnb,
                                                                         mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_sseb,
                                                                         mod4 = list_is_log_ev_ssir))

#SSIR
vec_post_probs_is_ssir = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssir,
                                                                         mod2 = list_is_log_ev_base, mod3 = list_is_log_ev_sseb,
                                                                         mod4 = list_is_log_ev_ssnb))
                                               

#PLOT
model_ev_method = 'IS'
data_type = 'NZ Waitemata 08/21 Subset I'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSEB = vec_post_probs_is_sseb, BASE = vec_post_probs_is_base,
                                                 SSNB = vec_post_probs_is_ssnb, SSI = vec_post_probs_is_ssir),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')

#************************
# FOUR MODELS #2
#************************
#BASE
vec_post_probs_is_base2 = GET_AGG_POSTERIOR_PROB(FLAG_BASELINE = TRUE, list_log_mod_evid = list(mod1 = list_is_log_ev_base2,
                                                                                               mod2 = list_is_log_ev_sseb2, mod3 = list_is_log_ev_ssnb2,
                                                                                               mod4 = list_is_log_ev_ssir2))

#SSEB
vec_post_probs_is_sseb2 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_sseb2,
                                                                         mod2 = list_is_log_ev_base2, mod3 = list_is_log_ev_ssnb2,
                                                                         mod4 = list_is_log_ev_ssir2))

#SSNB
vec_post_probs_is_ssnb2 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssnb2,
                                                                         mod2 = list_is_log_ev_base2, mod3 = list_is_log_ev_sseb2,
                                                                         mod4 = list_is_log_ev_ssir2))

#SSIR
vec_post_probs_is_ssir2 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_is_log_ev_ssir2,
                                                                         mod2 = list_is_log_ev_base2, mod3 = list_is_log_ev_sseb2,
                                                                         mod4 = list_is_log_ev_ssnb2))

#PLOT
model_ev_method = 'IS'
data_type = 'NZ Waitemata 08/21 II'
par(mfrow = c(2,1))
BOX_PLOT_POSTERIOR_PROBS(list_vec_results = list(SSEB = vec_post_probs_is_sseb2, BASE = vec_post_probs_is_base2,
                                                 SSNB = vec_post_probs_is_ssnb2, SSI = vec_post_probs_is_ssir2),
                         data_type = data_type, model_ev_method = '') #IS Model Evidence.')

#DATA
plot.ts(data_wait_08_21_sub2, ylab = 'Infection count',
        main = paste0('Waitemata, August 2021. Days ', t1, ' to ', t2))
