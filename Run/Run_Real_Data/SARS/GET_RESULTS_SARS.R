#**************************************
#* 
#* SARS RESULTS CANADA 2003
#* 
#***********************************
#RESULTS_FOLDER
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/2_SARS/RESULTS/WAVE_1/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/2_SARS/RESULTS/WAVE_2/'

#****************************************
#* GET MCMC FILES
#******************************************
list_results = GET_DATA_FILES(RESULTS_FOLDER)

#MCMC RESULTS
#R0
get_cis_and_mean(list_results$list_mcmc$mcmc_baseline$r0_vec)
get_cis_and_mean(list_results$list_mcmc$mcmc_sse$sse_params_matrix[,1])
get_cis_and_mean(list_results$list_mcmc$mcmc_ssi$ssi_params_matrix[,1])
get_cis_and_mean(list_results$list_mcmc$mcmc_sseb$r0_vec)
get_cis_and_mean(list_results$list_mcmc$mcmc_ssib$ssib_params_matrix[,1])

#MCMC EXTRACT
mcmc_baseline = list_results$list_mcmc$mcmc_baseline
mcmc_sse = list_results$list_mcmc$mcmc_sse
mcmc_ssi = list_results$list_mcmc$mcmc_ssi
mcmc_sseb = list_results$list_mcmc$mcmc_sseb
mcmc_ssib = list_results$list_mcmc$mcmc_ssib

#***************************************************************  
# 
# 2. MODEL COMPARISON RESULTS
#
#***************************************************************  
vec_mod_ev = GET_MODEL_EVIDENCE(list_results$list_mcmc, epidemic_data)

post_probs1 = GET_MODELS_POSTERIOR_PROBS(vec_mod_ev1, num_models = 5,
                                        prior_probs_models = c(0.2, 0.2, 0.2, 0.2, 0.2))

#SAVE
file_name = 'vec_mod_ev_sars1.rds'
saveRDS(vec_mod_ev1, file = paste0(RESULTS_FOLDER, file_name))
vec_mod_ev1 = readRDS(file = paste0(RESULTS_FOLDER, file_name))

file_name = 'vec_mod_ev_sars2.rds'
saveRDS(vec_mod_ev2, file = paste0(RESULTS_FOLDER, file_name))
vec_mod_ev2 = readRDS(file = paste0(RESULTS_FOLDER, file_name))

post_probs2 = GET_MODELS_POSTERIOR_PROBS(vec_mod_ev, num_models = 5,
                                        prior_probs_models = c(0.2, 0.2, 0.2, 0.2, 0.2))

#***************************************************************  
# 
# 4. MODEL PARAMS - BEST FITTING MODELS 
#
#***************************************************************  
#SSI Model; 
k = mcmc_ssi$ssi_params_matrix[,2]
#r0
get_cis_and_mean(mcmc_ssi$ssi_params_matrix[,1])
#SSI Model; k
get_cis_and_mean(mcmc_ssi$ssi_params_matrix[,2])

#SSIB Model
#r0
r0 = mcmc_ssib$ssib_params_matrix[,1]
a = mcmc_ssib$ssib_params_matrix[,2]
b = mcmc_ssib$ssib_params_matrix[,3]

#r0
get_cis_and_mean(mcmc_ssib$ssib_params_matrix[,1])
#a
get_cis_and_mean(mcmc_ssib$ssib_params_matrix[,2])
#b
get_cis_and_mean(mcmc_ssib$ssib_params_matrix[,3])

#PLOT
#HIST ALPHA
hist(k, breaks = 100, freq = FALSE,
     col = MODEL_COLORS[3],
     xlab = 'k', lwd = 1.5,
     cex.lab = 1.3, cex.main =1.3,
     main = bquote(paste('k, best fit SSI model, SARS Wave 2, 2003')))
