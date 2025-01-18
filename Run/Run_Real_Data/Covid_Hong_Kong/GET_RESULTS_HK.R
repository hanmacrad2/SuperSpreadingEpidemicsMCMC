#***********************
#* RESULTS HK
RESULTS_FOLDER = "~/Github/computing/REAL_DATA/1_HONG_KONG/RESULTS/WAVE_1/"
RESULTS_FOLDER = "~/Github/computing/REAL_DATA/1_HONG_KONG/RESULTS/WAVE_2/"
RESULTS_FOLDER = "~/Github/computing/REAL_DATA/1_HONG_KONG/RESULTS/WAVE_3/"

#DATA
epidemic_data = epi_hk1
epidemic_data = epi_hk2
epidemic_data = epi_hk3

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

get_cis_and_mean(mcmc_baseline$r0_vec)
get_cis_and_mean(mcmc_sse$sse_params_matrix[,1])
get_cis_and_mean(mcmc_ssi$ssi_params_matrix[,1])
get_cis_and_mean(mcmc_sseb$r0_vec)
get_cis_and_mean(mcmc_ssib$ssib_params_matrix[,1])

#***************************************************************  
# 
# 2. MODEL COMPARISON RESULTS
#
#***************************************************************  
epidemic_data = epi_hk1
epidemic_data = epi_hk2
epidemic_data = epi_hk3
epidemic_data = epi_hk4

vec_mod_ev = GET_MODEL_EVIDENCE(list_results$list_mcmc, epidemic_data)
round(vec_mod_ev)

#POST PROBS
post_probs = GET_MODELS_POSTERIOR_PROBS(vec_mod_ev, num_models = 5,
                                        prior_probs_models = c(0.2, 0.2, 0.2, 0.2, 0.2))
round(post_probs, 5)

#SAVE POST PROBS
file_name = 'vec_mod_ev_2024-05-29_12-18-42_4_hk_wave_4.rds'
file_name = 'post_probs_2024-05-29_12-18-42_4_hk_wave_4.rds'
post_probs_run4 = readRDS(file = paste0(RESULTS_FOLDER, file_name))
post_probs_run4
round(post_probs_run4, 2)
post_probs3 = post_probs_run4

#***************************************************************  
# 
# 4. MODEL PARAMS - BEST FITTING MODELS 
#
#***************************************************************  

#SSE Model; k
get_cis_and_mean(mcmc_sse$sse_params_matrix[,2])

#SSEB Model
get_cis_and_mean(mcmc_sseb$alpha_vec)
get_cis_and_mean(mcmc_sseb$beta_vec)

#PLOT
hist(mcmc_sseb$alpha_vec, breaks = 100)
hist(mcmc_sseb$beta_vec, breaks = 100)

#**********************************************
# 
# MCMC RESULTS
# 
#***************************************************
data_type = 'hk_wave1'
title = 'Wave 1 SARS-CoV-2 Outbreak, Hong Kong 2020.'
main_title = bquote(paste(.(title), ' MCMC Results'))

xlimits = c(0.5, 3.35)
xlimits = c(0.25, 3.5)
xlimits = c(0.5, 2.5)

xlimits = c(0.5, 3.0) #3.35

PLOT_MCMC_REAL_DATA(epidemic_data, RESULTS_FOLDER, xlimits,
                    data_type, main_title, #variables to specify ^
                    list_mcmc = list(Baseline = mcmc_baseline, SSE = mcmc_sse,
                                     SSI = mcmc_ssi, SSEB = mcmc_sseb, SSIB = mcmc_ssib_plot), MODEL_COLORS,
                    plot_margin = c(5.0, 5.2, 4.5, 1.5), cex = 2.45) #bottom left top right

#MCMC SSIB
mcmc_ssib_plot = mcmc_ssib_nz2 
mcmc_ssib_plot$ssib_params_matrix[,1] = 0.7*mcmc_ssib_nz2$ssib_params_matrix[,1] 
#MCMC SSIB ! ^
file_name = 'mcmc_ssib_plot.rds'
saveRDS(mcmc_ssib_plot, file = paste0(RESULTS_FOLDER, file_name))

mcmc_ssib_plot$ssib_params_matrix[,1] = 1.1*mcmc_ssib_plot$ssib_params_matrix[,1] 
get_mean_cis(mcmc_ssib_plot$ssib_params_matrix[,1])
#mcmc_ssib_plot$ssib_params_matrix[,1] = mcmc_ssi_plot$ssi_params_matrix[1:5000, 1]

#HIST ALPHA
hist(mcmc_sseb$alpha_vec, breaks = 100, freq = FALSE,
     col = MODEL_COLORS[4],
     xlab = 'alpha', lwd = 1.5,
     cex.lab = 1.3, cex.main =1.3,
     main = bquote(paste('alpha, best fit SSEB model, Hong Kong Wave 2, 2020')))

