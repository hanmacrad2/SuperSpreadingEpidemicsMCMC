#**************************#****************************************************
#*
#GET POSTERIOR PREDICTIVE PLOTS - HK
#
#******************************************************************************

#DATA
epidemic_data = epi_hk2

epidemic_data = df_nz_waita$cases
plot.ts(epidemic_data)

data_type = '  SARS-CoV-2 Outbreak, Wave 2 Hong Kong 2020'
DATA_SET = 'SARS-CoV-2, Wave 2 Hong Kong'

#*******************************************************************************
#* GET MCMC FILES
#******************************************v************************************
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/1_HONG_KONG/RESULTS/WAVE_2/'

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


#*********************************************************************************************
#*
#2. GET POSTEROR PREDICTIONS OF ALL MODELS
#
#*********************************************************************************************

#Posterior predictions
matrix_sim_baseline = SAMPLE_BASELINE_MCMC(mcmc_baseline, epidemic_data)
matrix_sim_sse = SAMPLE_SSE_MCMC(mcmc_sse, epidemic_data)
matrix_sim_ssi = SAMPLE_SSI_MCMC(mcmc_ssi, epidemic_data)
matrix_sim_sseb = SAMPLE_SSEB_MCMC(mcmc_sseb, epidemic_data)
matrix_sim_ssib = SAMPLE_SSIB_MCMC(mcmc_ssib, epidemic_data)

#********************
#PLOT
#********************

#SSE FIRST
FLAG_MODEL = GET_FLAGS_MODELS(SSE = TRUE)
title = paste0(' Models fit to ', DATA_SET)
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_sse, epidemic_data, df_hk_wave2,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[2],
                           RESULTS_FOLDER, plot_width = 11.0, plot_height = 8.0)

#GET_CI_PRED_DATA
GET_CI_PRED_DATA(matrix_sim_baseline, df_hk_wave2, MODEL_COLORS[1])
GET_CI_PRED_DATA(matrix_sim_ssi, df_hk_wave2, MODEL_COLORS[3]) #,CI = list(CI_95 = TRUE, CI_99 = FALSE)) #,  upper_quant = 0.99)
GET_CI_PRED_DATA(matrix_sim_sseb, df_hk_wave2, MODEL_COLORS[4])
GET_CI_PRED_DATA(matrix_sim_ssib, df_hk_wave2, MODEL_COLORS[5]) #, upper_quant = 0.66)
dev.off()

#GET_CI_PRED_DATA_SSIB(matrix_sim_ssi, df_hk_wave2, MODEL_COLORS[3], upper_quant = 0.99) 