#Apply model evidence - model comparison basic

#RUN MULTIPLE MCMC ITERATIONS

library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"

#***********************
#* METHOD 1: HARMONIC MEAN
#**********************

#***********************
# 1. RUN BASE MCMC
#**********************
list_log_ev_base = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER, FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE,
                                                                               SSNB = FALSE, SSIB = FALSE, SSIC = FALSE))
#Plot RESULTS
PLOT_MODEL_EV_RESULTS(list_log_ev_base, model_type = 'Baseline')
mean(list_log_ev_base)
sd(list_log_ev_base)

#***********************
# 2. RUN SSEB MCMC
#**********************
list_log_ev_sseb = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER) 
mean(list_log_ev_sseb)
sd(list_log_ev_sseb)

#Plot
PLOT_MODEL_EV_RESULTS(list_log_ev_sseb) 

#***********************
# 3. RUN SSNB MCMC
#**********************
list_log_ev_ssnb = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER, FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE,
                                                                               SSNB = TRUE, SSIB = FALSE, SSIC = FALSE))
mean(list_log_ev_ssnb)
sd(list_log_ev_ssnb)

#Plot
PLOT_MODEL_EV_RESULTS(list_log_ev_ssnb)

#***********************
# 4. RUN SSIB MCMC (Not correct)
#**********************
list_log_ev_ssib = LOAD_MCMC_GET_MODEL_EV_HM(OUTER_FOLDER, FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE,
                                                                               SSNB = FALSE, SSIB = TRUE, SSIC = FALSE))
mean(list_log_ev_ssib)
sd(list_log_ev_ssib)

#***********************
# 5. RUN SSIB MCMC
#**********************

#Plot
PLOT_MODEL_EV_RESULTS(list_log_ev_ssib)

#***********************
# 4. BAYES FACTORS
#**********************
log_bf_hm = list_log_ev_base - list_log_ev_sseb

#PLOT
PLOT_BAYES_FACTORS(log_bf_hm)
mean(log_bf_hm)

#***********************
# 5. POSTERIOR MODEL PROBABILITIES
#**********************
par(mfrow = c(1,1))
#GET_AGGREGATE_POSTERIOR_MODEL_PROB

#MODEL 1: BASE
vec_post_probs_hm = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_log_ev_base,
                                                            mod2 = list_log_ev_sseb, mod3 = list_log_ev_ssnb))

mean(vec_post_probs_hm); sd(vec_post_probs_hm)
PLOT_BAYES_FACTORS(vec_post_probs_hm)

#MODEL 2: SSEB
vec_post_probs_hm2 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_log_ev_sseb,
                                                                    mod2 = list_log_ev_base, mod3 = list_log_ev_ssnb))
mean(vec_post_probs_hm2); sd(vec_post_probs_hm2)
PLOT_BAYES_FACTORS(vec_post_probs_hm2)

#MODEL 3: SSNB
vec_post_probs_hm3 = GET_AGG_POSTERIOR_PROB(list_log_mod_evid = list(mod1 = list_log_ev_ssnb,
                                                                     mod2 = list_log_ev_base, mod3 = list_log_ev_sseb))
mean(vec_post_probs_hm3); sd(vec_post_probs_hm3)
PLOT_BAYES_FACTORS(vec_post_probs_hm3)

#MODEL 4: SSIB


#DATAFRAME OF POSTERIOR RESULTS
list_pp_results = list(Baseline_model = vec_post_probs_hm, SSEB_model = vec_post_probs_hm2, SSNB_model = vec_post_probs_hm3)
df_pp <- as.data.frame(do.call(cbind, list_pp_results))
df_pp
boxplot(df_pp, main = 'Posterior Model Probabilities (Model evidence via Harmonic Mean). Data - Baseline Model',
        col = c('red', 'green', 'blue'),
        cex.lab=1.3, cex.axis=1.3, cex.main=1.2, cex.sub=1.3)

#HISTOGRAM (OVERLAPPING?) OF RESULTS
hist(df_pp, main = 'Posterior Model Probabilities (Model evidence via Harmonic Mean). Data - Baseline Model')


#***********************
#* COMPARE WITH METHOD 2: IMPORTANCE SAMPLING -- MODEL EVIDENCE RESULTS
#**********************

#1. BASE
model_type = 'BASE'; run_number = 1
FOLDERX = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
lme_base_is = readRDS(file = paste0(FOLDERX, 'phat_ests_base_', run_number))

#Plot
mean(lme_base_is)
sd(lme_base_is)
PLOT_MODEL_EV_RESULTS(lme_base_is)

#2. SSEB
#INSPECT MODEL EVIDENCE RESULTS
model_type = 'SSEB'; run_number = 1
FOLDERX = paste0(OUTER_FOLDER, model_type, '/run_', run_number, '/')
lme_sseb_is = readRDS(file = paste0(FOLDERX, 'phat_ests_sseb_', run_number))

#Plot
mean(lme_sseb_is)
sd(lme_sseb_is)
PLOT_MODEL_EV_RESULTS(lme_sseb_is)

#PLOT TOTAL RESULTS
BOX_PLOT_MODEL_EV_RESULTS(list_log_ev_base, lme_base_is, list_log_ev_sseb, lme_sseb_is)

num_models = 3
for(i in 1:(num_models-1)){
  print(i)
}