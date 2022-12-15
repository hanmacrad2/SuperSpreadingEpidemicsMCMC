#PLOT MCMC GRID
library(coda)

PLOT_SSEC_MCMC_GRID <- function(epidemic_data, mcmc_output, n_mcmc,
                                seed_count, log_like_sim,
                                simulated = list(m1 = 1.2, m2 = 0.16),
                                mcmc_specs = list(model_type = 'Simulated',
                                                  mod_start_points = list(m1 = 1.2, m2 = 0.16), 
                                                  mod_par_names = c('R0', 'k'),
                                                  burn_in_pc = 0.05, thinning_factor = 10),
                                priors_list = list(r0_prior = c(1,0), k_prior = c(1,0)),
                                FLAGS_LIST = list(BURN_IN = FALSE, THIN = TRUE, PRIOR = TRUE,
                                                  ADAPTIVE = TRUE, MULTI_ALG = TRUE, PLOT_ADAPTIVE = FALSE)){
  
  #PLOT
  plot.new()
  par(mfrow=c(4,3))
  
  #EXTRACT MCMC SAMPLES
  log_like_mcmc = mcmc_output$log_like_vec; log_like_mcmc = unlist(log_like_mcmc)
  
  if (FLAGS_LIST$MULTI_ALG){
    m1_mcmc = mcmc_output$ssec_params_matrix[,1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output$ssec_params_matrix[,2]; m2_mcmc = unlist(m2_mcmc); m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
    
  } else {
    m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc);  m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
  }
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_specs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    mcmc_vec_size = mcmc_specs$n_mcmc
    print(paste0('mcmc vec size = ', mcmc_vec_size))
  }
  
  #BURN IN
  if (FLAGS_LIST$BURN_IN){
    burn_in = mcmc_specs$burn_in_pc*mcmc_vec_size
    m1_mcmc = m1_mcmc[burn_in:mcmc_vec_size]; m2_mcmc = m2_mcmc[burn_in:mcmc_vec_size]
    log_like_mcmc = log_like_mcmc[burn_in:mcmc_vec_size]
  } else {
    burn_in = 0
  }
  
  #LIMITS
  m1_min =  min(mcmc_specs$mod_start_points$m1, min(m1_mcmc, na.rm = TRUE));  m1_max =  max(mcmc_specs$mod_start_points$m1, max(m1_mcmc, na.rm = TRUE))
  m2_min = min(mcmc_specs$mod_start_points$m2, min(m2_mcmc, na.rm = TRUE)); m2_max = max(mcmc_specs$mod_start_points$m2, max(m2_mcmc, na.rm = TRUE))
  minll = min(min(log_like_mcmc, na.rm = TRUE), log_like_sim); maxll = max(max(log_like_mcmc,  na.rm = TRUE), log_like_sim)
  #LIMITS
  m1_lim = c(m1_min, m1_max);  m2_lim = c(m2_min, m2_max); lim_ll = c(minll, maxll)
  #Priors
  m2_prior =  paste0('exp(', priors_list$k_prior[1], ')')
  # m1_prior =  paste0('exp(', priors_list$a_prior_exp[1], ')')
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #************************
  #ROW 1: MCMC TRACE PLOTS
  #************************
  
  #r0
  r0_title = bquote(bold(R[0] ~ "MCMC. Simulated: " ~ .(simulated$m1)))
  plot.ts(m1_mcmc, ylab = mcmc_specs$mod_par_names[1], ylim= m1_lim, #bquote("Hello" ~ r[xy] == .(cor) ~ "and" ~ B^2)
          main = r0_title,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = simulated$m1, col = 'red', lwd = 2) #True = green
  
  if (FLAGS_LIST$PLOT_ADAPTIVE){
    sig1 = mcmc_output$sigma$sigma1_vec
    plot.ts(m1_mcmc, ylab = paste0(mcmc_specs$mod_par_names[1], ",sigma"), #ylim=c(min(min(sig1),min(m1_mcmc)), max(m1_mcmc)),
            main = paste(mcmc_specs$mod_par_names[1], "MCMC",
                         "Start: ", simulated$m1, ', Sigma (red)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma1_vec, col = 'red')
  }
  
  #***************
  #k
  plot.ts(m2_mcmc, ylab = mcmc_specs$mod_par_names[3], ylim= m2_lim,
          main = paste(mcmc_specs$mod_par_names[2], "MCMC",
                       "Simulated: ", simulated$m2),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = mcmc_specs$mod_start_points$m2, col = 'blue', lwd = 2) #True = green
  if (FLAGS_LIST$PLOT_ADAPTIVE){
    plot.ts(m2_mcmc, ylab = paste0(mcmc_specs$mod_par_names[2], ",sigma"), #ylim=c(0, m2_lim),
            main = paste(mcmc_specs$mod_par_names[2], "MCMC",
                         "Start: ", simulated$m2, ', Sigma (blue)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma2_vec, col = 'blue')
  }
  
  #***************
  #LOG LIKELIHOOD
  plot.ts(log_like_mcmc, ylab = "log likelihood", ylim= lim_ll,
          main = paste("Log Likelihood. N MCMC:", n_mcmc, ". Burn-in:", burn_in),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = log_like_sim, col = 'orange', lwd = 2)
  
  #**********************************************************
  #ROW 2:  HISTOGRAMS OF PARARMS (r0, k, loglike)
  #************************************************************
  
  #***********
  #HIST m1
  r0_titleII = bquote(bold(R[0] ~ "Prior: exp("~.(priors_list$r0_prior[1])~")"))
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_specs$mod_par_names[1], #ylab = 'Density',
       main = r0_titleII,
       xlim = m1_lim,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_specs$mod_start_points$m1, col = 'red', lwd = 2)
  
  #PRIOR PLOT 
  if (FLAGS_LIST$PRIOR) {
    xseq = seq(0, 1.5, length.out = 500)
    lines(xseq, dexp(xseq, priors_list$r0_prior[1]),
          type = 'l', lwd = 2, col = 'red')
  } else {
    #m2_prior = paste0('exp(', priors_list$b_prior_ga[1], ')')
  }
  
  #***********
  #HIST m2
  ktitle = bquote(bold(k ~ "Prior: exp("~.(priors_list$k_prior[1])~")"))
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_specs$mod_par_names[2], #ylab = 'Density',
       main = ktitle,
       xlim= m2_lim,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_specs$mod_start_points$m2, col = 'blue', lwd = 2)
  
  #Prior plot
  if (FLAGS_LIST$PRIOR) {
    xseq = seq(0, 1.5, length.out = 500)
    lines(xseq, dexp(xseq, priors_list$k_prior[1]),
          type = 'l', lwd = 2, col = 'blue')
  }
  
  #PRIOR PLOT
  # if (FLAGS_LIST$B_PRIOR_GAMMA) {
  #   xseq = seq(0, 0.3, length.out = 500)
  #   lines(xseq, dgamma(xseq, shape =  priors_list$b_prior_ga[1], scale =  priors_list$b_prior_ga[2]),
  #         type = 'l', lwd = 2, col = 'blue')
  # } else {
  #   xseq = seq(0, 10, length.out = 5000)
  #   lines(xseq, dexp(xseq, priors_list$b_prior_exp[1]),
  #         type = 'l', lwd = 2, col = 'blue')
  # }
  
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
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc)
  plot(seq_along(m1_mean), m1_mean,
       ylim= m1_lim,
       xlab = 'Time', ylab =  mcmc_specs$mod_par_names[1],
       main = bquote(bold(R[0] ~ "MCMC mean, Start:" ~ .(mcmc_specs$mod_start_points$m1))),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = mcmc_specs$mod_start_points$m1, col = 'red', lwd = 2)
  
  #m2 mean
  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  plot(seq_along(m2_mean), m2_mean,
       ylim= m2_lim,
       xlab = 'Time', ylab = mcmc_specs$mod_par_names[2],
       main = paste(mcmc_specs$mod_par_names[2], "MCMC mean, Start:", mcmc_specs$mod_start_points$m2),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  abline(h = mcmc_specs$mod_start_points$m2, col = 'blue', lwd = 2)
  
  #SCALING VEC
  plot.ts(mcmc_output$scaling_vec,  ylab = paste0('adaptive scaling vec'), 
          main = paste0('Adaptive scaling vec'),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #*****************
  #ROW 4: DATA INFECTIONS + JOINT DISTRIBUTIONS/MARGINALS
  #********************
  
  #i. TOTAL INFECTIONS
  inf_tite = paste0(mcmc_specs$model_type, " Data. Seed = ", seed_count) 
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = inf_tite,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #r0 VS K
  plot(m1_mcmc, m2_mcmc,
       xlab = mcmc_specs$mod_par_names[1], ylab = mcmc_specs$mod_par_names[2],
       main = paste0(mcmc_specs$mod_par_names[1], ' vs ', mcmc_specs$mod_par_names[2]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)
  
  #********************
  #v. DATAFRAME: RESULTS
  #********************
  
  df_results <- data.frame(
    rep = seed_count,
    n_mcmc = n_mcmc,
    mcmc_vec_size = mcmc_vec_size,
    r0_sim = simulated$m1[[1]],
    r0_start = mcmc_specs$mod_start_points$m1[[1]],
    r0_mean_mcmc = round(mean(m1_mcmc), 2), #round(mean(m1_mcmc[(mcmc_vec_size/2): mcmc_vec_size]), 2),
    k_sim = simulated$m2[[1]], 
    k_start =  mcmc_specs$mod_start_points$m2[[1]],
    k_mean_mcmc = round(mean(m2_mcmc), 2),
    log_like_sim = round(log_like_sim, 2),
    accept_rate = round(mcmc_output$accept_rate, 2),
    r0_es = round(effectiveSize(as.mcmc(m1_mcmc))[[1]], 2),
    k_es = round(effectiveSize(as.mcmc(m2_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) #format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
  
  print(df_results)
  
  return(df_results)
  
}
