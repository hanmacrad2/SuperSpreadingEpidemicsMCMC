#SET SARS P VALUES
#DATA
DATA_FOLDER = "~/Github/computing/REAL_DATA/2_SARS/DATA/"

#WAVE 1
data_type = '  SARS incidence data, Canada, Wave 1 2003'
DATA_SET = 'SARS, Canada Wave 1 2020'
file_name = 'df_sars_wave1.rds'
df_sars_wave1 = readRDS(file = paste0(DATA_FOLDER, file_name))

#EPI DATA
epidemic_data = df_sars_wave1$cases
plot.ts(epidemic_data)

#WAVE 2
data_type = '  SARS incidence data, Canada, Wave 2 2003'
DATA_SET = 'SARS, Canada Wave 2 2020'
file_name = 'df_sars_wave2.rds'
df_sars_wave2 = readRDS(file = paste0(DATA_FOLDER, file_name))

#EPI DATA
epidemic_data = df_sars_wave2$cases
plot.ts(epidemic_data)

#*******************************************************************************
#* GET MCMC FILES
#******************************************v************************************
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/2_SARS/RESULTS/WAVE_1/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/2_SARS/RESULTS/WAVE_2/'

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

#******************************************************************************
#
# GET P VALUES - MODEL CRITICISM
#
#******************************************************************************

#1. BASELINE
FLAGS_MODELS <- list(BASELINE = TRUE, SSE = FALSE, SSI = FALSE,
                     SSEB = FALSE, SSIB = FALSE)

p_values_baseline = RUN_MODEL_CRITICISM(epidemic_data, mcmc_baseline, FLAGS_MODELS = FLAGS_MODELS)

#2. SSE
FLAGS_MODELS <- list(BASELINE = FALSE, SSE = TRUE, SSI = FALSE,
                     SSEB = FALSE, SSIB = FALSE)

p_values_sse = RUN_MODEL_CRITICISM(epidemic_data, mcmc_sse, FLAGS_MODELS = FLAGS_MODELS)

#****************
#3. SSI
FLAGS_MODELS <- list(BASELINE = FALSE, SSE = FALSE, SSI = TRUE,
                     SSEB = FALSE, SSIB = FALSE)

p_values_ssi = RUN_MODEL_CRITICISM(epidemic_data, mcmc_ssi, FLAGS_MODELS = FLAGS_MODELS)

#4. SSEB 
FLAGS_MODELS <- list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                     SSEB = TRUE, SSIB = FALSE)
p_values_sseb = RUN_MODEL_CRITICISM(epidemic_data, mcmc_sseb, FLAGS_MODELS = FLAGS_MODELS)

#5. SSIB 
FLAGS_MODELS <- list(BASELINE = FALSE, SSE = FALSE, SSI = FALSE,
                     SSEB = FALSE, SSIB = TRUE)
p_values_ssib = RUN_MODEL_CRITICISM(epidemic_data, mcmc_ssib, FLAGS_MODELS = FLAGS_MODELS)


#CALL
p_values_baseline
p_values_sse
p_values_ssi
p_values_sseb
p_values_ssib

#*************
#SAVE LIST OF P VALUES
#P VALUES SARS 1
file_name = 'p_values_sars1.rds'
p_values_sars1 = list(p_values_baseline = p_values_baseline, p_values_sse = p_values_sse, 
                    p_values_ssi = p_values_ssi, p_values_sseb = p_values_sseb, p_values_ssib = p_values_ssib)

saveRDS(p_values_sars1, file = paste0(RESULTS_FOLDER, file_name))
p_values_sars1 = readRDS(file = paste0(RESULTS_FOLDER, file_name))

p.adjust(p_values_sars1$p_values_baseline, method = c("holm"))
p.adjust(p_values_sars1$p_values_ssib, method = c("holm"))

#P VALUES SARS 2
file_name = 'p_values_sars2.rds'
p_values_sars2 = list(p_values_baseline = p_values_baseline, p_values_sse = p_values_sse, 
                    p_values_ssi = p_values_ssi, p_values_sseb = p_values_sseb, p_values_ssib = p_values_ssib)

saveRDS(p_values_sars2, file = paste0(RESULTS_FOLDER, file_name))

#******************
#* SAVE P VALUES
file_name = 'p_values_ssib_sars_wave1.rds'
saveRDS(p_values_ssib, file = paste0(RESULTS_FOLDER, file_name))

file_name = 'p_values_baseline_sars_wave1.rds'
saveRDS(p_values_baseline, file = paste0(RESULTS_FOLDER, file_name))

file_name = 'p_values_sse_sars_wave1.rds'
saveRDS(p_values_sse, file = paste0(RESULTS_FOLDER, file_name))

file_name = 'p_values_ssi_sars_wave1.rds'
saveRDS(p_values_ssi, file = paste0(RESULTS_FOLDER, file_name))

file_name = 'p_values_sseb_sars_wave1.rds'
saveRDS(p_values_sseb, file = paste0(RESULTS_FOLDER, file_name))