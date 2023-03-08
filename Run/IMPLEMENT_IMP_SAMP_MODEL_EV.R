#IMPLEMENT MODEL EVIDENCE VIA IMPORTANCE SAMPLING 
library(SuperSpreadingEpidemicsMCMC)
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"

#***********************
# 1. BASE DATA (RUN AUTOMATICALLY)
#**********************

BASE_DATA_LOC = paste0(OUTER_FOLDER, 'BASE_DATA/')
data_baseI = readRDS(file = paste0(BASE_DATA_LOC, 'epi_data_base_1.rds'))
plot.ts(data_baseI)
runX = 1

#***************************
# 2. LOAD MCMC & GET MULTIPLE PHAT (log)
#***************************
OUTPUT_FOLDER = paste0(BASE_DATA_LOC, 'SSEB/')
ests_phat_sseb = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_LIST = list(BASE = FALSE, SSEB = TRUE,
                                                                                  SSIB = FALSE, SSIC = FALSE))
#SAVE (ADDED TO FUNCTION)
#saveRDS(ests_phat_sseb, file = paste0(OUTPUT_FOLDER, '/run_', runX, '/ests_phat_sseb.rds'))
ests_phat_sseb = readRDS(file = paste0(OUTPUT_FOLDER, 'run_', runX, '/phat_ests_sseb_vec100.rds')) 

OUTPUT_FOLDER = paste0(BASE_DATA_LOC, 'BASE/')
ests_phat_base = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_LIST = list(BASE = TRUE, SSEB = FALSE,
                                                       SSIB = FALSE, SSIC = FALSE))

OUTPUT_FOLDER = paste0(BASE_DATA_LOC, 'SSIB/')
ests_phat_ssib = LOAD_MCMC_GET_P_HAT(data_baseI, OUTPUT_FOLDER,
                                     FLAGS_LIST = list(BASE = FALSE, SSEB = FALSE,
                                                       SSIB = TRUE, SSIC = FALSE))

#***************************
# 3. GET POSTERIOR MODEL PROBABILITIES (PLOT) USE MULTIPLE PHATS
#***************************
post_probs_base = GET_AGG_POSTERIOR_PROB(list_log_phats = list(mod1 = ests_phat_base, mod2 = ests_phat_sseb,
                                                               mod3 = ests_phat_ssib))

post_probs_sseb = GET_AGG_POSTERIOR_PROB(list_log_phats = list(mod1 = ests_phat_sseb, mod2 = ests_phat_base, mod3 = ests_phat_ssib))

post_probs_ssib = GET_AGG_POSTERIOR_PROB(list_log_phats = list(mod1 = ests_phat_ssib, mod2 = ests_phat_base,
                                                               mod3 = ests_phat_sseb))
#******************************
#* PLOT POSTERIOR MODEL RESULTS
#*****************************
par(mfrow = c(2,1))

#BASE (log)
PLOT_MODEL_EV_RESULTS(post_probs_base)
#BASE
PLOT_MODEL_EV_RESULTS(post_probs_base, FLAG_RESULT_TYPE = list(phat = FALSE, post_prob = TRUE,
                                                                log = FALSE))
#SSEB
PLOT_MODEL_EV_RESULTS(post_probs_sseb, model_type = 'SSEB')
PLOT_MODEL_EV_RESULTS(post_probs_sseb,  model_type = 'SSEB', FLAG_RESULT_TYPE = list(phat = FALSE, post_prob = TRUE,
                                                               log = FALSE))

#SSIB
PLOT_MODEL_EV_RESULTS(post_probs_ssib, model_type = 'SSIB')
PLOT_MODEL_EV_RESULTS(post_probs_ssib,  model_type = 'SSIB', FLAG_RESULT_TYPE = list(phat = FALSE, post_prob = TRUE,
                                                               log = FALSE))

#******************
#* PLOT PHAT LOG RESULTS
#*****************
par(mfrow = c(2,1))

boxplot(ests_phat_sseb,
        ylab = 'Phat estimate (log) for SSEB',
        main = 'Phat estimates (log) for SSEB model. Base data. 100 reps')

hist(ests_phat_sseb, breaks = 50, freq = FALSE,
     xlab = 'Phat estimate (log) for SSEB',
     main = 'Phat estimates (log) for SSEB model. Base data. 100 reps')

#Base
boxplot(ests_phat_base,
        ylab = 'Phat estimate (log) for Baseline',
        main = 'Phat estimates (log) for Baseline model. Base data. 100 reps')

hist(ests_phat_base, breaks = 50, freq = FALSE,
     xlab = 'Phat estimate (log) for Baseline',
     main = 'Phat estimates (log) for Baseline model. Baseline data. 100 reps')

#SSIB
boxplot(ests_phat_ssib,
        ylab = 'Phat estimate (log) for SSIB',
        main = 'Phat estimates (log) for SSIB model. Base data. 100 reps')

hist(ests_phat_ssib, breaks = 50, freq = FALSE,
     xlab = 'Phat estimate (log) for SSIB',
     main = 'Phat estimates (log) for SSIB model. Base data. 100 reps')
