#GET P VALUES - NZ 

#DATA
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'

#DATA 1 NZ SOUTH, 2020
file_name = 'df_nz_south_20.rds'
df_nz_south = readRDS(file= paste0(DATA_FOLDER, file_name))
plot(df_nz_south$dates, df_nz_south$cases, type = 'l')

epidemic_data = df_nz_south$cases

#******************************
#DATA 2 NZ - WAITEMATA, AUCKLAND, 2021
file_name = 'df_nz_waita_21.rds'
df_nz_waita = readRDS(file= paste0(DATA_FOLDER, file_name))
plot(df_nz_waita$dates, df_nz_waita$cases, type = 'l')

epidemic_data = df_nz_waita$cases

#*******************************************************************************
#* GET MCMC FILES
#******************************************v************************************
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/NZ_WEDDING/run_2/'

RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/NZ_WAITEMATA/run_2/'

list_results = GET_DATA_FILES(RESULTS_FOLDER)

#MCMC RESULTS
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

#************************************************************************
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

#SAVE LIST OF P VALUES
#P VALUES NZ 1
file_name = 'p_values_nz1.rds'
p_values_nz1 = list(p_values_baseline = p_values_baseline, p_values_sse = p_values_sse, 
                   p_values_ssi = p_values_ssi, p_values_sseb = p_values_sseb, p_values_ssib = p_values_ssib)

saveRDS(p_values_nz1, file = paste0(RESULTS_FOLDER, file_name))
p_values_nz1 = readRDS(file = paste0(RESULTS_FOLDER, file_name))

#P VALUES NZ 2
file_name = 'p_values_nz2.rds'
p_values_nz2 = list(p_values_baseline = p_values_baseline, p_values_sse = p_values_sse, 
                    p_values_ssi = p_values_ssi, p_values_sseb = p_values_sseb, p_values_ssib = p_values_ssib)

saveRDS(p_values_nz2, file = paste0(RESULTS_FOLDER, file_name))
p_values_nz2 = readRDS(file = paste0(RESULTS_FOLDER, file_name))


#P VALUES
p_values_sse = c(0.08, 0.03, 0.05, 0.13, 0.64, 0.14, 0.17, 0.80, 0.21, 0.10)
p.adjust(p_values_sse, method = c("holm"))

#****************
#SAVE P VALUES
file_name = 'p_values_baseline_nz_1.rds'
saveRDS(p_values_baseline, file = paste0(RESULTS_FOLDER, file_name))

file_name = 'p_values_sse_nz_1.rds'
saveRDS(p_values_sse, file = paste0(RESULTS_FOLDER, file_name))

#SSI
file_name = 'p_values_ssi_nz_1.rds'
saveRDS(p_values_ssi, file = paste0(RESULTS_FOLDER, file_name))

#SSEB
file_name = 'p_values_sseb_nz_1.rds'
saveRDS(p_values_sseb, file = paste0(RESULTS_FOLDER, file_name))

#SSIB
file_name = 'p_values_ssib_nz_1.rds'
saveRDS(p_values_ssib, file = paste0(RESULTS_FOLDER, file_name))