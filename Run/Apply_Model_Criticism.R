#APPLY MODEL CRITICISM
library(zoo)
library(devtools)
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

#REP
rep = 1
print(paste0('rep = ', rep))
FOLDER_REP = paste0(RESULTS_FOLDER, "/rep_", rep, '/')
print(paste0('FOLDER_REP = ', FOLDER_REP))

#***********
#MCMC RESULTS
mcmc_output <- readRDS(paste0(FOLDER_REP, 'mcmc_output_rep_', rep, '.rds'))

#BASE MODEL
R0 = mcmc_output1$r0_vec[10]
print(R0)
sim_data = SIMULATE_BASELINE_EPIDEMIC(R0, num_days = 50)

#RESULTS
list_p_vals <- readRDS(paste0(results_rep, 'list_p_vals_rep', rep, '.rds'))

df_p_vals =  readRDS(paste0(RESULTS_FOLDER, '/df_total_p_values.rds'))


#Plot p values
PLOT_P_VALUES(df_p_vals, 'Baseline')

#OUTPUT 
list_p_vals <- readRDS(paste0(FOLDER_REP, 'list_p_vals_rep', rep, '.rds'))
df_p_vals <- readRDS(paste0(FOLDER_REP, 'df_total_p_values.rds'))
df_p_vals
df_summary_stats <- readRDS(paste0(FOLDER_REP, 'df_replicated_summary_stats_', rep, '.rds'))
df_summary_stats
df_true_sum_stats <- readRDS(paste0(FOLDER_REP, 'df_ss_true_rep_', rep, '.rds')) 
df_true_sum_stats
list_ss_iters_i1 <- readRDS(paste0(FOLDER_REP, 'list_ss_iters_i1.rds'))
list_ss_iters_i1
length(list_ss_iters_i1)
len_data = length(list_p_vals)
len_data

#DISPLAY SUMMARY STATS
PLOT_SUMMARY_STATS(FOLDER_REP, canadaX, data_type, rep)

PLOT_REPLICATED_DATA(FOLDER_REP, canadaX, rep, data_type)
