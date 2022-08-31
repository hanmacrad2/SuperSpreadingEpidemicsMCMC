#APPLY MODEL CRITICISM
library(zoo)
library(SuperSpreadingEpidemicsMCMC)
#library(coda)

#FOLDER
ROOT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism_II/"
ITER = 'iter_1/'
MODEL_TYPE = 'baseline'
RESULTS_FOLDER =  paste0(ROOT_FOLDER, ITER, MODEL_TYPE)
print(RESULTS_FOLDER)

#MCMC
modelling_specs = list(n_reps = 10, n_mcmc = 10000)

#APPLY
RUN_MODEL_CRITICISM(canadaX, RESULTS_FOLDER, modelling_specs = modelling_specs)

#Run
rep = 1
mcmc_output1 <- readRDS(paste0(RESULTS_FOLDER, '/rep_1/mcmc_output_rep_', rep, '.rds'))

#REP
print(paste0('rep = ', rep))
folder_rep = paste0(root_folder, "/rep_", rep, '/')
print(paste0('folder_rep = ', folder_rep))

#MCMC
mcmc_output <- readRDS(paste0(folder_rep, 'mcmc_output_rep_', rep, '.rds'))

#BASE MODEL
R0 = mcmc_output1$r0_vec[10]
print(R0)
sim_data = SIMULATE_BASELINE_EPIDEMIC(R0, num_days = 50)
