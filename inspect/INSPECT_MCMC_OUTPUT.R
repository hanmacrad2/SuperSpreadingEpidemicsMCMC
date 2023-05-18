#INSPECT MCMC OUTPUT

#SSNB
run = 5
model_type = 'ssnb'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')

#MCMC
i = 1
mcmc_output_nb_2 = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
#mcmc_samples =  cbind(mcmc_output$ssic_params_matrix, mcmc_output$eta_matrix)

#PLOT
plot.ts(mcmc_output_nb_2$ssnb_params_matrix)

#SSIR
run = 5
model_type = 'ssir'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')

#MCMC
i = 1
mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
#mcmc_samples =  cbind(mcmc_output$ssic_params_matrix, mcmc_output$eta_matrix)

#PLOT
plot.ts(mcmc_output$ssir_params_matrix)
plot.ts(mcmc_output$eta_matrix)

#COMPARE
library(EpiEstim)
help(estimate_R)
