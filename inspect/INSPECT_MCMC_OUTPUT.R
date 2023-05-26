#INSPECT MCMC OUTPUT

#FOLDERS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/"

OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSNB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIR_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSEB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIB_DATA/')

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

#SSIB
run = 1
model_type = 'ssib'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')

#MCMC
i = 1
mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
#mcmc_samples =  cbind(mcmc_output$ssic_params_matrix, mcmc_output$eta_matrix)

#PLOT
plot.ts(mcmc_output$a_vec)
plot.ts(mcmc_output$b_vec)
plot.ts(mcmc_output$c_vec)
plot.ts(mcmc_output$r0_vec)
plot.ts(mcmc_output$log_like_vec)
plot.ts(mcmc_output$sigma$sigma5)
plot(mcmc_output$non_ss)


plot.ts(mcmc_output$eta_matrix)


#COMPARE
library(EpiEstim)
help(estimate_R)
