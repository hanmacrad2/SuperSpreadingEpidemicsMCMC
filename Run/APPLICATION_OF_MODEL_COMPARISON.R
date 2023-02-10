#APPLY MODEL COMPARISON

#*************************************
#* PART II:  SSEB Model vs BASE MODEL
#************************************

#DATA
data_sseb1 = SIMULATE_EPI_SSEB()
plot.ts(data_sseb1)

#SAVE DATA
saveRDS(data_sseb1, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb', seedX, '.rds' ))

#MODEL EVIDENCE
list_model_ev_sseb1 = RUN_MODEL_EV_SSEB(data_sseb1, CURRENT_OUTPUT_FOLDER)
saveRDS(list_model_ev_sseb1, file = paste0(CURRENT_OUTPUT_FOLDER, '/list_model_ev_sseb.rds' ))

#PLOT LOG EVIDENCE
par(mfrow = c(2,1))
plot(seq_along(list_model_ev_sseb1), list_model_ev_sseb1,
     ylim = c(min(list_model_ev_sseb1)-50, max(list_model_ev_sseb1)+50), lwd = 1, pch = 19,
     main = 'SSEB Model - Model Evidence(log). SSEB data',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')

#PLOT 2 (Close up)
plot(seq_along(list_model_ev_sseb1), list_model_ev_sseb1,
     ylim = c(min(list_model_ev_sseb1)-0.5, max(list_model_ev_sseb1)+0.5), lwd = 1, pch = 19,
     main = 'SSEB Model - Model Evidence(log) (Close-up). SSEB data',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')

#***************************
#* LOG EVIDENCE -- SSEB DATA -- BASELINE MODEL
#***************************

#MODEL EVIDENCE
list_model_ev_base2 = RUN_MODEL_EV_BASE(data_sseb1, CURRENT_OUTPUT_FOLDER)
saveRDS(list_model_ev_base2, file = paste0(CURRENT_OUTPUT_FOLDER, '/list_model_ev_base2.rds' ))

#PLOT LOG EVIDENCE
par(mfrow = c(1,1))
plot(seq_along(list_model_ev_base2), list_model_ev_base2,
     ylim = c(min(list_model_ev_base2)-50, max(list_model_ev_base2)+50), lwd = 1, pch = 19,
     main = 'Baseline - Model Evidence(log). SSEB data',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')

#PLOT 2 (Close up)
plot(seq_along(list_model_ev_base2), list_model_ev_base2,
     ylim = c(min(list_model_ev_base2)-0.5, max(list_model_ev_base2)+0.5), lwd = 1, pch = 19,
     main = 'Baseline Model - Model Evidence(log) (Close-up). SSEB data',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')


#***************************
#* BAYES FACTORS *TO DO
#***************************
log_bfs2 = list_model_ev_sseb1 - list_model_ev_base2 #LOG_MODEL_EVIDENCE(mcmc_sse_output$log_like_vec)
bayes_factors2 = exp(log_bfs2)
bayes_factors2_total = bayes_factors2

bayes_factors2 = bayes_factors2_total[2:30]

#PLOT BFS
hist(bayes_factors2, breaks = 200,
     lwd = 1, pch = 19, freq = FALSE,
     main = 'Bayes Factors. SSEB vs Baseline. SSEB data',
     xlab = 'Bayes Factor')
axis(side=1, at=seq(0,200, 5), labels=seq(0,200,5))

plot(seq_along(bayes_factors2), bayes_factors2,
     lwd = 1, pch = 19, 
     ylim = c(min(bayes_factors2)-2, max(bayes_factors2)+2), 
     main = 'Bayes Factors. SSEB vs Baseline. SSEB. data',
     xlab = 'mcmc iteration', ylab = 'Bayes Factor')

#***************************
#* RJMCMC 
#***************************

#RUN RJMCMC
n_mcmc = 50000
start_time = Sys.time()
print(paste0('start_time:', start_time))
rj_sseb1 = RJMCMC_BASE_SSEB(data_sseb1, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
#SAVE
saveRDS(rj_sseb1, file = paste0(CURRENT_OUTPUT_FOLDER, '/rj_sseb1', seedX, '.rds' ))

#PLOT 
true_r0 = 1.8
df_sseb = PLOT_SSEB_RJMCMC(data_sseb1, rj_sseb1, n_mcmc)

#RUN MULTIPLE

#*****************
#SSEB
#*****************
#SSEB CREATE OUTPUT FOLDER

OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/rjmcmc_rjmcmc_ssebbase1"

CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

rj_m2 = RUN_RJMCMC_MULT(data_sseb1, CURRENT_OUTPUT_FOLDER)

#*****************
#BASE
#*****************
#BASE CREATE OUTPUT FOLDER
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/rjmcmc_base1"

CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)


rj_m1 = RUN_RJMCMC_MULT(dataI, CURRENT_OUTPUT_FOLDER)
