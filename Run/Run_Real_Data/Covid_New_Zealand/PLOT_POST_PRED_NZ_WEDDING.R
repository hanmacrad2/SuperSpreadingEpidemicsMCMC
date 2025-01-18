#**************************#****************************************************
#*
#GET POSTERIOR PREDICTIVE PLOTS
#
#******************************************************************************
#DATA
file_name = 'df_nz.rds'
df_nz = readRDS(file= paste0(DATA_FOLDER, file_name))

#EPIDEMIC DATA
epidemic_data = df_nz$cases
plot.ts(epidemic_data)

data_type = '  SARS-CoV-2 - New Zealand South Outbreak, 2020'
DATA_SET = 'SARS-CoV-2, NZ South, 2020'


#*******************************************************************************
#* GET MCMC FILES
#******************************************v************************************
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/NZ_WEDDING/run_2/'

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

#*******************************
#GET POSTERIOR PREDICTIONS

#*********************************************************************************************
#*
#1. PLOT POSTEROR PREDICTIONS OF ALL MODELS
#
#*********************************************************************************************

#Posterior predictions
matrix_sim_baseline = SAMPLE_BASELINE_MCMC(mcmc_baseline, epidemic_data)
matrix_sim_sse = SAMPLE_SSE_MCMC(mcmc_sse, epidemic_data)
matrix_sim_ssi = SAMPLE_SSI_MCMC(mcmc_ssi, epidemic_data) #, SSI = TRUE)
matrix_sim_sseb = SAMPLE_SSEB_MCMC(mcmc_sseb, epidemic_data)
matrix_sim_ssib = SAMPLE_SSIB_MCMC(mcmc_ssib, epidemic_data)

#********************
#PLOT
#********************

#SSIB FIRST
FLAG_MODEL = GET_FLAGS_MODELS(SSI = TRUE)
title = paste0(' Models fit to ', DATA_SET)
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_ssi, epidemic_data, df_nz,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[3],
                           RESULTS_FOLDER, EPI_DATA = FALSE, ALL_MODELS = TRUE,
                           plot_width = 11.0, plot_height = 7.5, cex = 1.25) #, upper_quant = 0.99)

#GET_CI_PRED_DATA
GET_CI_PRED_DATA(matrix_sim_baseline, df_nz, MODEL_COLORS[1])# , upper_quant = 0.99)
GET_CI_PRED_DATA(matrix_sim_sse, df_nz, MODEL_COLORS[2]) #,  upper_quant = 0.99)
#GET_CI_PRED_DATA(matrix_sim_ssi, df_nz, MODEL_COLORS[3], upper_quant = 0.99)
GET_CI_PRED_DATA(matrix_sim_sseb, df_nz, MODEL_COLORS[4]) #,  upper_quant = 0.99)
GET_CI_PRED_DATA(matrix_sim_ssib, df_nz, MODEL_COLORS[5]) #,  upper_quant = 0.69)

dev.off()

#CI = list(CI_95 = FALSE, CI_99 = TRUE)
#SSI FIRST
FLAG_MODEL = GET_FLAGS_MODELS(SSI = TRUE)
title = paste0(' Models fit to ', DATA_SET)
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_ssi, epidemic_data, df_sars_wave2,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[3],
                           RESULTS_FOLDER, EPI_DATA = FALSE, ALL_MODELS = TRUE,
                           plot_width = 11.0, plot_height = 7.5, cex = 1.25)

#GET_CI_PRED_DATA
GET_CI_PRED_DATA(matrix_sim_baseline, df_sars_wave2, MODEL_COLORS[1])
GET_CI_PRED_DATA(matrix_sim_sse, df_sars_wave2, MODEL_COLORS[2])
GET_CI_PRED_DATA(matrix_sim_sseb, df_sars_wave2, MODEL_COLORS[4])
GET_CI_PRED_DATA(matrix_sim_ssib, df_sars_wave2, MODEL_COLORS[5])
dev.off()

#**********************
#* MODELS SEPARATELY
#BASELINE 
matrix_sim_baseline = SAMPLE_BASELINE_MCMC(mcmc_baseline, epidemic_data)

#PLOT
FLAGS_MODELS = GET_FLAGS_MODELS(BASELINE = TRUE)
title = paste0('BASELINE model fit to', DATA_SET)
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_baseline, epidemic_data, df_sars_wave2,
                           title, data_type, FLAGS_MODELS, MODEL_COLORS[1],
                           RESULTS_FOLDER, EPI_DATA = FALSE)

#SSE 
matrix_sim_sse = SAMPLE_SSE_MCMC(mcmc_sse, epidemic_data)

#PLOT
FLAG_MODEL = GET_FLAGS_MODELS(SSE = TRUE)
title = 'SSE model fit to SARS-CoV-2, HK Wave 1 2020' 
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_sse, epidemic_data, df_hk_wave1,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[2],
                           RESULTS_FOLDER, EPI_DATA = FALSE)

#SSI
matrix_sim_ssi = SAMPLE_SSE_MCMC(mcmc_ssi, epidemic_data, SSI = TRUE)

#PLOT
FLAG_MODEL = GET_FLAGS_MODELS(SSI = TRUE)
title = 'SSI model fit to SARS-CoV-2, HK Wave 1 2020'
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_ssi, epidemic_data, df_hk_wave1,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[3],
                           RESULTS_FOLDER, EPI_DATA = FALSE)

#*************
#SSEB
matrix_sim_sseb = SAMPLE_SSEB_MCMC(mcmc_sseb, epidemic_data)

#PLOT
FLAG_MODEL = GET_FLAGS_MODELS(SSEB = TRUE)
title = 'SSEB model fit to SARS-CoV-2, HK Wave 1 2020'
data_type = ' SARS-CoV-2 incidence data, Hong Kong Wave 1 2020'
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_sseb, epidemic_data, df_hk_wave1,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[4],
                           RESULTS_FOLDER, EPI_DATA = FALSE)

#*************
#SSIB
matrix_sim_ssib = SAMPLE_SSIB_MCMC(mcmc_ssib, epidemic_data)

#PLOT
FLAG_MODEL = GET_FLAGS_MODELS(SSIB = TRUE)
title = 'SSIB model fit to SARS-CoV-2, HK Wave 1 2020'
POSTERIOR_PREDICTIVE_PLOTS(matrix_sim_ssib, epidemic_data, df_hk_wave1,
                           title, data_type, FLAG_MODEL, MODEL_COLORS[5],
                           RESULTS_FOLDER, EPI_DATA = FALSE)
