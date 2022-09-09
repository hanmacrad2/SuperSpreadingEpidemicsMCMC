#APPLY MODEL CRITICISM
library(zoo)
library(devtools)
library(SuperSpreadingEpidemicsMCMC)
#library(coda)

#FOLDER
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism_II/"
ITER = 'iter_22/'
MODEL_TYPE = 'baseline'; DATA_TYPE = 'Canadian'
ROOT_FOLDER =  paste0(OUTER_FOLDER, ITER, MODEL_TYPE, '/')
print(ROOT_FOLDER)
#MCMC
modelling_specs = list(n_reps = 1000, n_mcmc = 50000)

#*************************************************
#* 1. MODEL CRITICISM
#APPLY MCMC
RUN_MODEL_CRITICISM(canadaX, ROOT_FOLDER, modelling_specs = modelling_specs)

#**************************************************
#* 2. DISPLAY MODEL CRITICISM RESULTS
DISPLAY_MODEL_CRITICISM(ROOT_FOLDER, canadaX, MODEL_TYPE, DATA_TYPE, REP)

#****************************************************
#DISPLAY RESULTS

#TOTAL P VALUES 
PLOT_P_VALUES(ROOT_FOLDER, MODEL_TYPE)

#REPLICATED DATA (SPECIFIC RESULTS)
REP = 500; print(paste0('REP = ', REP))
FOLDER_REP = paste0(ROOT_FOLDER, "/REP_", REP, '/')
#SUMMARY STATS
PLOT_SUMMARY_STATS(FOLDER_REP, canadaX, DATA_TYPE, REP)
#REPLICATED DATA
PLOT_REPLICATED_DATA(FOLDER_REP, canadaX, REP, DATA_TYPE)


#MCMC OUTPUT
mcmc_output <- readRDS(paste0(FOLDER_REP, 'mcmc_output_rep_', REP, '.rds'))
PLOT_BASELINE_R0_MCMC(canadaX, mcmc_output, DATA_TYPE)

#*********************************************************************
#OUTPUT - INSPECT
#*********************************************************************
FOLDER_REP = paste0(ROOT_FOLDER, "/REP_", REP, '/')
print(paste0('FOLDER_REP = ', FOLDER_REP))

#DISPLAY RESULTS
PLOT_P_VALUES(ROOT_FOLDER, MODEL_TYPE)
PLOT_SUMMARY_STATS(FOLDER_REP, canadaX, DATA_TYPE, REP)
PLOT_REPLICATED_DATA(FOLDER_REP, canadaX, REP, DATA_TYPE)

#OUTPUT 
list_p_vals <- readRDS(paste0(FOLDER_REP, 'list_p_vals_REP', REP, '.rds'))
df_p_vals <- readRDS(paste0(ROOT_FOLDER, 'df_total_p_values.rds'))
df_p_vals
df_summary_stats <- readRDS(paste0(FOLDER_REP, 'df_replicated_summary_stats_', REP, '.rds'))
df_summary_stats
df_true_sum_stats <- readRDS(paste0(FOLDER_REP, 'df_ss_true_rep_', REP, '.rds')) 
df_true_sum_stats
list_ss_iters_i1 <- readRDS(paste0(FOLDER_REP, 'list_ss_iters_i1.rds'))
list_ss_iters_i1
length(list_ss_iters_i1)
len_data = length(list_p_vals)
len_data

#*********************************************************************
#BRAINSTORM
#*********************************************************************
#P VALUES
df_p_values =  readRDS(paste0(ROOT_FOLDER, 'df_total_p_values.rds'))

#MCMC RESULTS
mcmc_output <- readRDS(paste0(FOLDER_REP, 'mcmc_output_rep_', REP, '.rds'))

#MCMC OUTPUT
PLOT_BASELINE_R0_MCMC(canadaX, mcmc_output, DATA_TYPE)

#BASE MODEL
R0 = mcmc_output1$r0_vec[10]
print(R0)
sim_data = SIMULATE_BASELINE_EPIDEMIC(R0, num_days = 50)

#RESULTS
list_p_vals <- readRDS(paste0(results_REP, 'list_p_vals_REP', REP, '.rds'))
list_p_vals
df_p_vals =  readRDS(paste0(ROOT_FOLDER, '/df_total_p_values.rds'))
PLOT_P_VALUES(df_p_vals, 'Baseline')

#*********************************************************************
# P VALUES
#*********************************************************************

#FOLDER REP
print(paste0('rep = ', REP))
FOLDER_REP = paste0(ROOT_FOLDER, "rep_", REP, '/')
print(paste0('FOLDER_REP = ', FOLDER_REP))

#GET TRUE SUMMARY STATISTICS
df_true_ss = readRDS(paste0(FOLDER_REP, 'df_ss_true_rep_', REP, '.rds' )) #RENAME!!
print('passed I')
#df_true_ss = get_summary_stats(true_rep_sim, TRUE)

#GET REPLICATED SUMMARY STATISTICS 
df_summary_stats_rep <- readRDS(paste0(FOLDER_REP, 'df_replicated_summary_stats_', REP, '.rds' ))
print('passed II')

#GET P VALUES
list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) GET_P_VALUE(df_true_ss[,x], df_summary_stats_rep[,x]))
list_p_vals

saveRDS(list_p_vals, file = paste0(FOLDER_REP, '/list_p_vals_rep', rep, ".rds"))

#SUMMARY STATS #2
mid_point = length(canadaX)/2
mid_point
end_point  = length(canadaX)
end_point
sum(canadaX[(mid_point+1):end_point])

sum(canadaX[1:mid_point])

