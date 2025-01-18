#GET P VALUES HK

#DATA
epidemic_data = epi_hk1

#*******************************************************************************
#* GET MCMC FILES
#******************************************v************************************
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/1_HONG_KONG/RESULTS/WAVE_1/'

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
p_values_hk = list(p_values_baseline = p_values_baseline, p_values_sse = p_values_sse, 
                   p_values_ssi = p_values_ssi, p_values_sseb = p_values_sseb, p_values_ssib = p_values_ssib)
file_name = 'p_values_hk1.rds'
saveRDS(p_values_hk, file = paste0(RESULTS_FOLDER, file_name))

p_values_hk = readRDS(file = paste0(RESULTS_FOLDER, file_name))

#LIST p values
p_values_hk$p_values_baseline
p_values_hk$p_values_sse
p_values_hk$p_values_ssi
p_values_hk$p_values_sseb
p_values_hk$p_values_ssib

p_values_hk_base <- unname(unlist(p_values_hk$p_values_baseline))
p_values_hk_sse <- unname(unlist(p_values_hk$p_values_sse))
p_values_hk_ssi <- unname(unlist(p_values_hk$p_values_ssi))
p_values_hk_sseb <- unname(unlist(p_values_hk$p_values_sseb))
p_values_hk_ssib <- unname(unlist(p_values_hk$p_values_ssib))

p.adjust(p_values_hk_sseb, method = c("holm"))

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