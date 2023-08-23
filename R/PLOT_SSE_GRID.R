#PLOT SSE MCMC GRID
PLOT_SSE_MCMC_GRID <- function(epidemic_data, mcmc_output, n_mcmc, 
                                sim_vals = list(m1 =0.16, m2 = 1.8),
                                mcmc_specs = list(model_type = 'Simulated',
                                                  mod_par_names = c('R0', 'k'),
                                                  burn_in_pc = 0.2, thinning_factor = 10),
                                FLAGS_LIST = list(PRIOR = TRUE, MULTI_ALG = TRUE, 
                                                  PLOT_ADAPTIVE = FALSE)){
  
  #PLOT
  plot.new(); par(mfrow=c(4,3))
  
  #PRIORS
  PRIORS_USED =  GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSE() 
  
  #PRIOR R0
  if(PRIORS_USED$SSE$R0$EXP){
    mr0_prior = paste0('exp(', list_priors$r0[1], ')')
    xseq_r0 = seq(0, 3, length.out = 500)
    dr0e = dexp(xseq_r0, list_priors$r0[1])
  }
  
  #PRIOR k 
  if(PRIORS_USED$SSE$k$EXP){
    m2_prior = paste0('exp(', list_priors$k[1], ')')
    x2 = seq(0, 3, length.out = 500)
    d2 = dexp(x2, list_priors$k[1])
  }
  
  #MCMC + LIKELIHOOD SAMPLES EXTRACT 
  sse_sim_params = c(sim_vals$m1, sim_vals$m2)
  lambda_vec =  get_lambda(epidemic_data); 
  log_like_sim = LOG_LIKE_SSNB(epidemic_data, lambda_vec, sse_sim_params) 
  log_like_mcmc = mcmc_output$log_like_vec; log_like_mcmc = unlist(log_like_mcmc)
  
  if (FLAGS_LIST$MULTI_ALG){
    m1_mcmc = mcmc_output$ssnb_params_matrix[,1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output$ssnb_params_matrix[,2]; m2_mcmc = unlist(m2_mcmc); m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
    
  } else {
    m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc);  m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
  }

  #LIMITS
  m1_min =  min(sim_vals$m1, min(m1_mcmc, na.rm = TRUE));  m1_max =  max(sim_vals$m1, max(m1_mcmc, na.rm = TRUE))
  m2_min = min(sim_vals$m2, min(m2_mcmc, na.rm = TRUE)); m2_max = max(sim_vals$m2, max(m2_mcmc, na.rm = TRUE))
  minll = min(min(log_like_mcmc, na.rm = TRUE), log_like_sim); maxll = max(max(log_like_mcmc,  na.rm = TRUE), log_like_sim)
  m1_lim = c(m1_min, m1_max);  m2_lim = c(m2_min, m2_max); lim_ll = c(minll, maxll)
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #*****************
  #ROW 1: DATA INFECTIONS + JOINT DISTRIBUTIONS/MARGINALS
  #********************
  
  #i. TOTAL INFECTIONS
  inf_tite = paste0(mcmc_specs$model_type, " Data. R0: ", sim_vals$m1, 'k:', sim_vals$m2) 
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = inf_tite,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #r0 VS K
  plot(m1_mcmc, m2_mcmc,
       xlab = mcmc_specs$mod_par_names[1], ylab = mcmc_specs$mod_par_names[2],
       main = paste0(mcmc_specs$mod_par_names[1], ' vs ', mcmc_specs$mod_par_names[2]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #Empty
  plot.ts(0, xlab = '', ylab = '')
  
  #************************
  #ROW 2: MCMC TRACE PLOTS
  #************************
  
  #R0
  r0_title = bquote(bold(R[0] ~ "MCMC. sim: " ~ .(sim_vals$m1)))
  plot.ts(m1_mcmc, ylab = mcmc_specs$mod_par_names[1], ylim = m1_lim, #bquote("Hello" ~ r[xy] == .(cor) ~ "and" ~ B^2)
          main = r0_title,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m1, col = 'red', lwd = 2) #True = green
  
  if (FLAGS_LIST$PLOT_ADAPTIVE){
    sig1 = mcmc_output$sigma$sigma1_vec
    plot.ts(m1_mcmc, ylab = paste0(mcmc_specs$mod_par_names[1], ",sigma"), #ylim=c(min(min(sig1),min(m1_mcmc)), max(m1_mcmc)),
            main = paste(mcmc_specs$mod_par_names[1], "MCMC",
                         "Start: ", sim_vals$m1, ', Sigma (red)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma1_vec, col = 'red')
  }
  
  #***************
  #k
  plot.ts(m2_mcmc, ylab = mcmc_specs$mod_par_names[2], ylim= m2_lim,
          main = paste(mcmc_specs$mod_par_names[2], "MCMC.",
                       "sim: ", sim_vals$m2),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m2, col = 'blue', lwd = 2) #True = green
  
  if (FLAGS_LIST$PLOT_ADAPTIVE){
    plot.ts(m2_mcmc, ylab = paste0(mcmc_specs$mod_par_names[2], ",sigma"), #ylim=c(0, m2_lim),
            main = paste(mcmc_specs$mod_par_names[2], "MCMC",
                         "Start: ", sim_vals$m2, ', Sigma (blue)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma2_vec, col = 'blue')
  }
  
  #***************
  #LOG LIKELIHOOD
  plot.ts(log_like_mcmc, ylab = "log likelihood", ylim= lim_ll,
          main = paste("Log Likelihood. N MCMC:", n_mcmc, ". Burn-in:", mcmc_specs$burn_in_pc),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = log_like_sim, col = 'orange', lwd = 2) 
  
  #**********************************************************
  #ROW 2:  HISTOGRAMS OF PARARMS (r0, k, loglike)
  #************************************************************
  
  #***********
  #HIST m1
  r0_titleII = bquote(bold(R[0])) # ~ "Prior: exp("~.(priors_list$r0_prior[1])~")"))
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_specs$mod_par_names[1], #ylab = 'Density',
       main = paste("R0, prior:", mr0_prior),
       xlim = m1_lim,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = sim_vals$m1, col = 'red', lwd = 2)
  
  #PRIOR PLOT 
  if(FLAGS_LIST$PRIOR) {
    lines(xseq_r0, dr0e, type = 'l', lwd = 2) 
  }
 
  #***********
  #HIST m2
  ktitle = bquote(bold(k)) # ~ "Prior: exp("~.(priors_list$k_prior[1])~")"))
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_specs$mod_par_names[2], #ylab = 'Density',
       main = paste(mod_par_names[2],
                    " prior:", m2_prior),
       xlim= m2_lim,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = sim_vals$m2, col = 'blue', lwd = 2)
  
  #PRIOR PLOT 
  if(FLAGS_LIST$PRIOR){
    lines(x2, d2, type = 'l', lwd = 2) 
  }
  
  #***********
  #HIST LOG_LIKE_VEC
  hist(log_like_mcmc, freq = FALSE, breaks = 100, xlim= lim_ll,
       xlab = 'Log likelihood', #ylab = 'Density',
       main = 'Log likelihood',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = log_like_sim, col = 'orange', lwd = 2)
  
  #************************************************
  #ROW 3: CUMULATIVE MEAN PLOTS
  #************************************************
  
  #m1 mean
  titleX = bquote(R[0] ~ "MCMC mean, Start:" ~ .(sim_vals$m1))
  PLOT_CUM_MEAN_MCMC(m1_mcmc, titleX = titleX, ylabX =  mcmc_specs$mod_par_names[1],
                     ylim = c(m1_min, m1_max))
  abline(h = sim_vals$m1, col = 'blue', lwd = 2)
  
  #m2 mean
  titleX = paste(mcmc_specs$mod_par_names[2], "MCMC mean, Start:", sim_vals$m2)
  PLOT_CUM_MEAN_MCMC(m2_mcmc, title = titleX, ylabX =  mcmc_specs$mod_par_names[2],
                     ylim = c(m2_min, m2_max))
  abline(h = sim_vals$m2, col = 'red', lwd = 2)
  
  #SCALING VEC
  plot.ts(mcmc_output$scaling_vec,  ylab = paste0('adaptive scaling vec'), 
          main = paste0('Adaptive scaling vec'),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #********************
  #v. DATAFRAME: RESULTS
  #********************
  
  df_results <- data.frame(
    #rep = seed_count,
    n_mcmc = n_mcmc,
    #n_thin_samps = n_samples,
    r0_sim = sim_vals$m1[1],
    r0_mean_mcmc = round(mean(m1_mcmc), 2), 
    r0_lo_95_cred_int = round(get_lower_ci(m1_mcmc), 2), 
    r0_up_95_cred_int = round(get_upper_ci(m1_mcmc), 2), 
    k_sim = sim_vals$m1[2], 
    k_mean_mcmc = round(mean(m2_mcmc), 2),
    k_lo_95_cred_int = round(get_lower_ci(m2_mcmc), 2),
    k_up_95_cred_int = round(get_upper_ci(m2_mcmc), 2), 
    log_like_sim = round(log_like_sim, 2),
    accept_rate = round(mcmc_output$accept_rate, 2),
    r0_es = round(effectiveSize(as.mcmc(m1_mcmc))[[1]], 2),
    k_es = round(effectiveSize(as.mcmc(m2_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) #format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
  
  #print(df_results)
  
  return(df_results)
  
}
