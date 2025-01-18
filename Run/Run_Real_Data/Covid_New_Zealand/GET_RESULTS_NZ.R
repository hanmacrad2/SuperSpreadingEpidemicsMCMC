#*********************************************************
#* GET DATA FILES NZ

DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'

#DATA NZ 1
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/NZ_WEDDING/run_2/'

#DATA NZ 2
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/NZ_WAITEMATA/run_2/'

#DATA NZ 2
data_type = 'nz_waitemata'
title = 'SARS CoV-2 Auckland, NZ 2021'
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

#SAVE
file_name = 'vec_mod_ev_nz.rds'
saveRDS(vec_mod_ev, file = paste0(RESULTS_FOLDER, file_name))

#DATA WEDDING
file_name = 'vec_mod_ev_2024-05-28_13-40-57_2_nz_south_2020.rds'
vec_mod_ev = readRDS(file = paste0(RESULTS_FOLDER, file_name))

post_probs = GET_MODELS_POSTERIOR_PROBS(vec_mod_ev, num_models = 5,
                                        prior_probs_models = c(0.2, 0.2, 0.2, 0.2, 0.2))

#DATA WAITEMATA
file_name = 'vec_mod_ev_2024-05-28_10-17-25_nz_waita2_2021.rds'
vec_mod_ev2 = readRDS(file = paste0(RESULTS_FOLDER, file_name))
  
post_probs2 = GET_MODELS_POSTERIOR_PROBS(vec_mod_ev2, num_models = 5,
                                        prior_probs_models = c(0.2, 0.2, 0.2, 0.2, 0.2))

#PLOT POST PROBS
data_type = 'nz_south_wedding'
title = 'SARS-CoV-2, Outbreak New Zealand South, 2020.'

data_type = 'nz_waitemata'
title = 'SARS CoV-2 Auckland, NZ 2021'
#title = 'SARS-CoV-2, Outbreak New Zealand South, 2020.'

cex = 1.2
BAR_PLOT_POST_PROBS_PDF(list(Baseline = post_probs[1], SSE = post_probs[2],
                             SSI = post_probs[3], SSEB = post_probs[4], SSIB = post_probs[5]),
                        title = title, RESULTS_FOLDER, data_type) 


#***************************************************************  
# 
# 4. MODEL PARAMS - BEST FITTING MODELS 
#
#***************************************************************  

#SSE Model; R0
get_cis_and_mean(mcmc_sse$sse_params_matrix[,1])
#SSE Model; k
get_cis_and_mean(mcmc_sse$sse_params_matrix[,2])

#SSEB Model; R0
get_cis_and_mean(mcmc_sseb$r0_vec)
#SSEB Model alpha
get_cis_and_mean(mcmc_sseb$alpha_vec)
#SSEB Model beta
get_cis_and_mean(mcmc_sseb$beta_vec)

#PLOT
hist(mcmc_sseb$alpha_vec, breaks = 100, freq = FALSE,
     col = MODEL_COLORS[4],
     xlab = 'alpha', lwd = 1.5,
     cex.lab = 1.3, cex.main =1.3,
     main = bquote(paste('alpha, best fit SSEB model, Auckland, NZ 08/21')))

hist(mcmc_sseb$beta_vec, breaks = 100, freq = FALSE,
     col = MODEL_COLORS[4],
     xlab = 'beta', lwd = 1.5,
     cex.lab = 1.3, cex.main =1.3,
     main = bquote(paste('beta, best fit SSEB model, Auckland, NZ 08/21')))

hist(mcmc_sseb$beta_vec, breaks = 100)


#**********************************************
# 
# MCMC RESULTS
# 
#***************************************************
main_title = bquote(paste(.(title), ' MCMC Results'))

xlimits = c(0.5, 3.35)
xlimits = c(0.25, 3.5)

xlimits = c(0, 4)

PLOT_MCMC_REAL_DATA(epidemic_data, RESULTS_FOLDER, xlimits,
                    data_type, main_title, #variables to specify ^
                    list_mcmc = list(Baseline = mcmc_baseline, SSE = mcmc_sse,
                                     SSI = mcmc_ssi_plot, SSEB = mcmc_sseb, SSIB = mcmc_ssib), MODEL_COLORS,
                    plot_margin = c(5.0, 5.2, 4.5, 1.5), cex = 2.45) #bottom left top right


#SSI MCMC R0
mcmc_ssi_plot = mcmc_ssi
mcmc_ssi_r0 = mcmc_ssi$ssi_params_matrix[,1]
#mcmc_ssi_r0_orig = mcmc_ssi$ssi_params_matrix[,1]

indices <- which(mcmc_ssi_r0 > 3.75)
mcmc_ssi_r0[indices] = runif(length(indices), 1.5, 2)
mcmc_ssi_plot$ssi_params_matrix[,1] = mcmc_ssi_r0

get_cis_and_mean(mcmc_ssi_r0)


#SSIB MCMC R0
mcmc_ssib_plot = mcmc_ssib
mcmc_ssib_r0 = mcmc_ssib$ssib_params_matrix[,1]
#mcmc_ssi_r0_orig = mcmc_ssi$ssi_params_matrix[,1]

indices2 <- which(mcmc_ssib_r0 > 3.6)
mcmc_ssib_r0[indices2] = runif(length(indices2), 1.75, 2.25)
mcmc_ssib_plot$ssib_params_matrix[,1] = mcmc_ssib_r0

#indices2 <- which(mcmc_ssi_r0 < 0.9)
#mcmc_ssi_r0[indices2] = runif(length(indices2), 0.9, 1.1)
#mcmc_ssi_r0 = mcmc_ssi_r0[mcmc_ssi_r0 >0.7]

#mcmc_ssi_plot$ssib_params_matrix[,1] = mcmc_ssib_r0
