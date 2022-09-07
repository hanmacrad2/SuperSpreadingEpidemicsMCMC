#APPLY MODEL CRITICISM
library(zoo)
library(devtools)
library(SuperSpreadingEpidemicsMCMC)
#library(coda)

#FOLDER
ROOT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism_II/"
ITER = 'iter_2/'
MODEL_TYPE = 'baseline'; DATA_TYPE = 'Canadian'
RESULTS_FOLDER =  paste0(ROOT_FOLDER, ITER, MODEL_TYPE, '/')
print(RESULTS_FOLDER)
#MCMC
modelling_specs = list(n_reps = 1000, n_mcmc = 50000)

#*************
#* 1. MODEL CRITICISM
#APPLY MCMC
RUN_MODEL_CRITICISM(canadaX, RESULTS_FOLDER, modelling_specs = modelling_specs)

#REP
REP = 1; print(paste0('REP = ', REP))
FOLDER_REP = paste0(RESULTS_FOLDER, "/REP_", REP, '/')

#*************
#* 2. DISPLAY MODEL CRITICISM RESULTS
DISPLAY_MODEL_CRITICISM(RESULTS_FOLDER, canadaX, MODEL_TYPE, DATA_TYPE, REP)

#*********************************************************************
#BRAINSTORM
#*********************************************************************
FOLDER_REP = paste0(RESULTS_FOLDER, "/REP_", REP, '/')
print(paste0('FOLDER_REP = ', FOLDER_REP))

#DISPLAY RESULTS
PLOT_P_VALUES(RESULTS_FOLDER, MODEL_TYPE)
PLOT_SUMMARY_STATS(FOLDER_REP, canadaX, DATA_TYPE, REP)
PLOT_REPLICATED_DATA(FOLDER_REP, canadaX, REP, DATA_TYPE)
PLOT_BASELINE_R0_MCMC(canadaX, mcmc_output, DATA_TYPE)

#OUTPUT 
list_p_vals <- readRDS(paste0(FOLDER_REP, 'list_p_vals_REP', REP, '.rds'))
df_p_vals <- readRDS(paste0(FOLDER_REP, 'df_total_p_values.rds'))
df_p_vals
df_summary_stats <- readRDS(paste0(FOLDER_REP, 'df_REPlicated_summary_stats_', REP, '.rds'))
df_summary_stats
df_true_sum_stats <- readRDS(paste0(FOLDER_REP, 'df_ss_true_REP_', REP, '.rds')) 
df_true_sum_stats
list_ss_iters_i1 <- readRDS(paste0(FOLDER_REP, 'list_ss_iters_i1.rds'))
list_ss_iters_i1
length(list_ss_iters_i1)
len_data = length(list_p_vals)
len_data

#MCMC RESULTS
mcmc_output <- readRDS(paste0(FOLDER_REP, 'mcmc_output_REP_', REP, '.rds'))

#BASE MODEL
R0 = mcmc_output1$r0_vec[10]
print(R0)
sim_data = SIMULATE_BASELINE_EPIDEMIC(R0, num_days = 50)

#RESULTS
list_p_vals <- readRDS(paste0(results_REP, 'list_p_vals_REP', REP, '.rds'))
df_p_vals =  readRDS(paste0(RESULTS_FOLDER, '/df_total_p_values.rds'))
PLOT_P_VALUES(df_p_vals, 'Baseline')

