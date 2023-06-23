#INSPECT MCMC OUTPUT

#FOLDERS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/"

OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSNB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIR_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSEB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'SSIB_DATA/')
OUTER_FOLDER = paste0(OUTER_FOLDER, 'MOCK_DATA/')

#PARAMS
NMCMC = 30000

#BASELINE
OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASELINE_DATA/DATA_BASELINE_2/')
run = '3_gp'
run = 1

model_type = 'baseline'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/'); print(CURRENT_FOLDER)

#MCMC
i = 10
mcmc_output_baseline_exp1 = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))

mcmc_output_baseline_gamma = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))

#PLOT
PLOT_BASELINE_R0_MCMC(data_baseline, mcmc_output_baseline_exp1) #exp(1) - MUCH TIGHTER!!

PLOT_BASELINE_R0_MCMC(data_baseline, mcmc_output_baseline_gamma)

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
model_type = 'ssir'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')

#MCMC
i = 1
mcmc_output_ssir = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
#mcmc_samples =  cbind(mcmc_output$ssic_params_matrix, mcmc_output$eta_matrix)

#PLOT
plot.ts(mcmc_output$ssir_params_matrix)
plot.ts(mcmc_output$eta_matrix)
PLOT_SSIR_MCMC_GRID(data_ssib, mcmc_output_ssir, 0, 1,
                                0, NMCMC)

#SSIB #GET DATA
run = 1
model_type = 'ssib'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')

#run gamma
run = 'run_2_ga_prior'; 
OUTER_FOLDER = paste0(OUTER_FOLDER, 'DATA_SSIB_2/') #SSIB_DATA
CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
CURRENT_FOLDER

#MCMC
i = 10
mcmc_output = readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', model_type, '_', i ,'.rds'))
mcmc_samples =  cbind(mcmc_output$a_vec, mcmc_output$b_vec,  mcmc_output$c_vec)
plot.ts(mcmc_samples)

#PLOT
plot.ts(mcmc_output$a_vec)
plot.ts(mcmc_output$b_vec)
plot.ts(mcmc_output$c_vec)
plot.ts(mcmc_output$r0_vec)
plot.ts(mcmc_output$log_like_vec)
plot.ts(mcmc_output$sigma$sigma5)
plot(mcmc_output$non_ss)

#PLOT GRID
PLOT_SSB_MCMC_GRID(data_ssib2, mcmc_output, 
                   sim_vals = list(m1 = 0.6, m2 = 0.1, m3 = 10),
                   FLAGS_MODELS = list(SSEB = FALSE, SSIB = TRUE))

#COMPARE
library(EpiEstim)
help(estimate_R)
