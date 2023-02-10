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
#* BAYES FACTORS *TO DO
#***************************
log_bfs1 = list_model_ev_base[1:50] - list_model_ev_sseb1 #LOG_MODEL_EVIDENCE(mcmc_sse_output$log_like_vec)
bayes_factors1 = exp(log_bfs1)

#PLOT BFS
hist(bayes_factors1, breaks = 200,
     lwd = 1, pch = 19, freq = FALSE,
     main = 'Bayes Factors. Baseline vs SSEB. Baseline data',
     xlab = 'Bayes Factor')
axis(side=1, at=seq(0,200, 5), labels=seq(0,200,5))

plot(seq_along(bayes_factors1), bayes_factors1,
     lwd = 1, pch = 19, 
     ylim = c(min(bayes_factors1)-2, max(bayes_factors1)+2), 
     main = 'Bayes Factors. Baseline vs SSEB. Baseline data',
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

#PLOT (*Include in plotting function: time_elap, If SSEB ALL TRUE PARAMS!)
true_r0 = 1.8
PLOT_SSB_MCMC_REAL_DATA(data_sseb1, rj_sseb1, n_mcmc, true_r0)