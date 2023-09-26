#GET DF RESULTS

#DATAFRAME
df_sse1 = df_results
file2 = 'sse_mcmc_2023-08-30_19-40-19.rds'
file3 = 'sse_mcmc_2023-08-31_17-02-05.rds'
df_sse3 = readRDS(paste0(COMP_FOLDER, file3))
df_results = df_sse3


#SUBSET TOTAL INFECTIONS SMALL
df_results_tot2 = df_results[df_results$tot_infs > max_infs, ]
df_results = df_results_tot2
