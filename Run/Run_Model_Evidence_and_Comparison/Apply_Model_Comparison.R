#Apply Model Comparison
library(coda)
library(SuperSpreadingEpidemicsMCMC)
plot.new()
par(mfrow=c(1,1))

#FOLDER
n_mcmc = 50000
seedX = 0
set.seed(seedX)
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/1_base_sseb"

#*****************
#BASE
#*****************
#BASE CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

#SSEB CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_sseb', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)


#*****************
#SSEB FOLDER
#*****************
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/2_sseb_base"
#SSEB CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

#******************************
#* 1. PART I: BASE MODEL vs SSEB Model
#*****************************
true_r0 = 1.6

#1. SIMULATE BASE MODEL DATA
data_baseI = SIMULATE_BASELINE_EPIDEMIC(R0 = true_r0)
plot.ts(data_baseI)

#SAVE DATA
saveRDS(data_baseI, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_base', seedX, '.rds' ))

#2.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_base_1 = MCMC_INFER_BASELINE(data_baseI)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_base_1$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_base_1, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base', seedX, '.rds' ))

#PLOT
true_loglike = LOG_LIKE_BASELINE(data_baseI, true_r0)
PLOT_BASELINE_R0_MCMC(data_baseI, mcmc_base_1, true_r0, true_loglike, data_type = 'Baseline')

#MODEL EVIDENCE
log_mod_ev1 = LOG_MODEL_EVIDENCE(mcmc_base_1$log_like_vec)
list_log_ev = c(log_mod_ev1)

#REPEAT MODEL EVIDENCE FOR BASE MODEL
RUN_MODEL_EV_BASE <- function(data_base, n_reps = 100){
  
  #List of model evidences
  list_log_ev = c()
  
  #REPEAT FOR REPS
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_base = MCMC_INFER_BASELINE(data_base)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_base$time_elap = time_elap
    saveRDS(mcmc_base, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_base', i, '.rds' ))
    
    #MODEL EVIDENCE
    log_mod_ev = LOG_MODEL_EVIDENCE(mcmc_base$log_like_vec)
    list_log_ev = c(list_log_ev, log_mod_ev)
    print(log_mod_ev)
  
}

  return(list_log_ev)
  
}

#APPLY
n_reps = 100
list_model_ev_base = RUN_MODEL_EV_BASE(dataI)
saveRDS(list_model_ev_base, file = paste0(CURRENT_OUTPUT_FOLDER, '/list_model_ev_base.rds' ))

#PLOT MODEL EVIDENCE
par(mfrow = c(2,1))

#PLOT 1
plot(seq_along(list_model_ev_base), list_model_ev_base,
     ylim = c(min(list_model_ev_base)-50, max(list_model_ev_base)+50), lwd = 1, pch = 19,
     main = 'Base Model - Model Evidence(log)',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')

#PLOT 2 (Close up)
plot(seq_along(list_model_ev_base), list_model_ev_base,
     ylim = c(min(list_model_ev_base)-0.5, max(list_model_ev_base)+0.5), lwd = 1, pch = 19,
        main = 'Base Model - Model Evidence(log) (Close-up)',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')

#**********************************************************************************
#2. SSEB

#RUN MCMC
mcmc_sse_output = MCMC_INFER_SSEB(dataI, n_mcmc)

#2. SSE
df_sse = PLOT_SS_MCMC_GRID(dataI, mcmc_sse_output, n_mcmc,
                           mcmc_specs = list(model_type = 'SSE', n_mcmc = 50000, 
                                             mod_start_points = list(m1 = 0.8, m2 = 0.1, m3 = 10),
                                             mod_par_names = c('alpha', 'beta', 'gamma'),
                                             seed_count = 1, burn_in_pc = 0.05, thinning_factor = 10),
                           priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                              c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                           FLAGS_LIST = list(SSI = FALSE, BURN_IN = TRUE, THIN = TRUE,
                                             DATA_AUG = TRUE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                             PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))


#EDIT FUNCTION & RUN
RUN_MODEL_EV_SSEB <- function(data_base, CURRENT_OUTPUT_FOLDER, n_reps = 50){
  
  #INITIALISE
  list_log_ev = c()
  
  for(i in 1:n_reps){
    
    print(paste0('i =', i))
    #RUN MCMC
    start_time = Sys.time()
    print(paste0('start_time:', start_time))
    mcmc_sseb = MCMC_INFER_SSEB(data_base, n_mcmc)
    end_time = Sys.time()
    time_elap = get_time(start_time, end_time)
    mcmc_sseb$time_elap = time_elap
    saveRDS(mcmc_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb', i, '.rds' ))
    
    #MODEL EVIDENCE
    log_mod_ev = LOG_MODEL_EVIDENCE(mcmc_sseb$log_like_vec)
    list_log_ev = c(list_log_ev, log_mod_ev)
    print(log_mod_ev)
    
  }
  
  return(list_log_ev)
}

#APPLY
list_model_ev_base_sseb = RUN_MODEL_EV_SSEB(dataI)
saveRDS(list_model_ev_base_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/list_model_ev_sseb_base.rds' ))

#PLOT LOG EVIDENCE
plot(seq_along(list_model_ev_base_sseb), list_model_ev_base_sseb,
     ylim = c(min(list_model_ev_base_sseb)-50, max(list_model_ev_base_sseb)+50), lwd = 1, pch = 19,
     main = 'SSEB Model - Model Evidence(log). Baseline data',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')

#PLOT 2 (Close up)
plot(seq_along(list_model_ev_base_sseb), list_model_ev_base_sseb,
     ylim = c(min(list_model_ev_base_sseb)-0.5, max(list_model_ev_base_sseb)+0.5), lwd = 1, pch = 19,
     main = 'SSEB Model - Model Evidence(log) (Close-up). Baseline data',
     xlab = 'mcmc iteration', ylab = 'log(model evidence)')


#***************************
#* Bayes Factors
#***************************
log_bfs1 = list_model_ev_base[1:50] - list_model_ev_base_sseb #LOG_MODEL_EVIDENCE(mcmc_sse_output$log_like_vec)
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

#*************************
#* RJMCMC_BASE_SSEB
#**************************
rj_out1 = RJMCMC_BASE_SSEB(dataI, n_mcmc)
bf = rj_out1$bf
#Fix
rj_out1$list_accept_rates$accept_rate_b = 0
rj_out1$list_accept_rates$accept_rate_g = 0

#PLOT
PLOT_SSB_MCMC_REAL_DATA(dataI, rj_out1, n_mcmc, true_r0)

#RUN ii
rj_out2 = RJMCMC_BASE_SSEB(dataI, n_mcmc)

#PLOT
PLOT_SSB_MCMC_REAL_DATA(dataI, rj_out2, n_mcmc, true_r0)


#*************************************
#* 1. PART II: BASE MODEL vs SSEB Model
#* ***********************************

#DATA
data_sseb1 = SIMULATE_EPI_SSEB()
plot.ts(data_sseb1)


#*************************
#* 3. SSEB VS SSIB Model
#************************

#1. RUN SSI MODEL
mcmc_ssi_output = SSI_MCMC_ADAPTIVE(sim_data_canadaX1)

#2. RUN SSE MODEL
mcmc_sse_output = SSE_MCMC_ADAPTIVE(canadaX)

#3. GET HARMONIC MEANS
bf = get_bayes_factor_harmonic_means(mcmc_ssi_output$log_like_vec, mcmc_sse_output$log_like_vec)
print(paste0('bayes factor = ', bf))

#****************************
#PLOT OUTPUT

#1. SSI
df_ssi = PLOT_SS_MCMC_GRID(sim_data_canadaX1, mcmc_ssi_output) 

#2. SSE
df_sse = PLOT_SS_MCMC_GRID(canadaX, mcmc_sse_output,
                              mcmc_specs = list(model_type = 'SSE', n_mcmc = 500000, 
                                                mod_start_points = list(m1 = 0.8, m2 = 0.1, m3 = 10),
                                                mod_par_names = c('alpha', 'beta', 'gamma'),
                                                seed_count = 1, burn_in_pc = 0.05, thinning_factor = 10),
                              priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                                 c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                              FLAGS_LIST = list(SSI = FALSE, BURN_IN = TRUE, THIN = TRUE,
                                                DATA_AUG = TRUE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                                PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))

#II. SIGMA
PLOT_SIGMA_ADADPTIVE(mcmc_ssi_output, mcmc_specs)
