#' Grid Plot of MCMC Super-Spreading Results
#'Grid Plot of MCMC Results for Super-Spreading Epidemic Models
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param mcmc_output mcmc samples from mcmc sampler/algorithm
#' @param mcmc_specs A list of mcmc specifications
#' \itemize{
#'   \item \code{"model_type"} - Model type; Super Spreading Individuals \code{'SSI'} or Super Spreading Events \code{'SSE'}
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler (integer)
#'   \item \code{"mod_start_points"} - Model parameter starting points; where the mcmc algorithm started sampling from
#'   \item \code{"mod_par_names"} - Names of the model parameters, e.g \code{"a, b, c"}
#'   \item \code{"seed_count"} - Seed for data generation & mcmc iteration
#'   \item \code{"burn_in_pc"} - Proportion of mcmc samples to remove at the start as burn-in
#'   \item \code{"thinning_factor"}  - factor of total \code{"n_mcmc"} size of which samples are kept. Only if  \code{"FLAGS_LIST$THIN = TRUE"}, otherwise all samples are kept
#' }
#' @param FLAGS_LIST A list of boolean variables for switching on/off certain functionality
#' \itemize{
#'   \item \code{"THIN"}  - Return a thinned mcmc sample if TRUE, reduced by a factor of \code{"thinning_factor"}
#'   \item \code{"PRIOR"}  - Plot prior distributions in grid if TRUE
#' }
#' @export
#'
#'@author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
PLOT_SSB_MCMC_GRID <- function(epidemic_data, mcmc_output, n_mcmc = 50000, 
                               sim_vals = list(m1 = 0.6, m2 = 0.1, m3 = 10), #0.8, 0.2, 10
                                         mcmc_specs = list(seed_count = 1,  burn_in_pc = 0.2, thinning_factor = 10),
                                        FLAGS_MODELS = list(SSEB = FALSE, SSIB = FALSE),
                                         FLAGS_LIST = list(THIN = TRUE, PRIOR = TRUE,
                                                           PLOT_SSIB_DATA = FALSE)){

  #PLOT
  plot.new()
  par(mfrow=c(5,5))
  
  #PRIORS
  PRIORS_USED =  GET_PRIORS_USED() 
  
  #MODEL
  if (FLAGS_MODELS$SSEB){
    model_type = 'SSE-B'
    mod_par_names = c('alpha', 'beta', 'gamma')
    
  } else if (FLAGS_MODELS$SSIB){
    model_type = 'SSI-B'
    mod_par_names = c('a', 'b', 'c')
    #PRIORS
    list_priors = GET_PRIORS_SSIB() 
  }
  
  #PRIORS
  if(PRIORS_USED$SSIB$a$BETA){
    m1_prior = paste0('Beta(', list_priors$a[1], ', ', list_priors$a[2], ')')
    xseq_a = seq(0, 1.5, length.out = 500)
    dba = dbeta(xseq_a, list_priors$a[1], list_priors$a[2])
  }
  
  if(PRIORS_USED$SSIB$c$GAMMA){
    m3_prior = paste0('Gamma(', list_priors$c[1], ', ', list_priors$c[2], ')')
    xseq_c = seq(0, 10, length.out = 500)
    dgc = dgamma(xseq_c, list_priors$c[1], list_priors$c[2])
  }
  
  if(PRIORS_USED$SSIB$R0$EXP){
    mr0_prior = paste0('exp(', list_priors$r0[1], ')')
    xseq_r0 = seq(0, 3, length.out = 500)
    dr0e = dexp(xseq_r0, list_priors$r0[1])
  }
  
  #MCMC SAMPLES
    m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc);  m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
    m3_mcmc = mcmc_output[3]; m3_mcmc = unlist(m3_mcmc);  m3_mcmc = m3_mcmc[!is.na(m3_mcmc)]
    r0_mcmc = mcmc_output[4]; r0_mcmc = unlist(r0_mcmc); r0_mcmc = r0_mcmc[!is.na(r0_mcmc)]
    log_like_mcmc = mcmc_output$log_like_vec; log_like_mcmc = unlist(log_like_mcmc)
    r0_sim = sim_vals$m1 + (sim_vals$m2*sim_vals$m3)
  
  #THINNING FACTOR
  if(FLAGS_LIST$THIN){
    thinning_factor = mcmc_specs$thinning_factor
    mcmc_vec_size = n_mcmc/thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    thinning_factor = 1; mcmc_vec_size = n_mcmc
    print(paste0('mcmc vec size = ', mcmc_vec_size))
  }
  
  #CUMULATIVE MEANS + PARAM SAMPLE LIMITS
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc);  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  m3_mean = cumsum(m3_mcmc)/seq_along(m3_mcmc); r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  
  #DATA (STARTING DATA)
  if (FLAGS_LIST$PLOT_SSIB_DATA){
    non_ss_start = epidemic_data[[1]]; ss_start = epidemic_data[[2]]
    epidemic_data = non_ss_start + ss_start
  }
  
  #******************************************************************
  #* PLOTS *

  #************************
  #ROW 1: DATA INFECTIONS
  #************************
  print(epidemic_data)
  #i. TOTAL INFECTIONS
  #inf_title = paste0(model_type, " Data") #mcmc_specs$seed_count, ' '
  plot.ts(epidemic_data, ylab = 'Daily Infections count',
          main = 'SSIB Data')
         # cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) # xlab = 'Time',
  
  if (FLAGS_LIST$PLOT_SSIB_DATA){
  #**********************
  #ii DATA - NON SS
  #**********************
  plot.ts(non_ss_start, ylab = 'count', main = 'Non SS, Start (blk), Final (aq), Mean (or)')
  lines(mcmc_output$non_ss[mcmc_vec_size, ], col = 'aquamarine', lwd = 2) #Final
  lines(colMeans(mcmc_output$non_ss), col = '#FF3300', lwd = 2) #mean

  #iii. MEAN T1, T2
  plot.ts(colMeans(mcmc_output$non_ss[((mcmc_vec_size/2)+1):mcmc_vec_size, ]), ylab = 'count',
          main = 'Non SS - Mean T1 (prp), Mean T2 (or)', col = 'orange')
  lines(colMeans(mcmc_output$non_ss[1:(mcmc_vec_size/2), ]), col = 'blueviolet', lwd = 2)

  #**********************
  #v, iv. DATA SS
  #**********************
  plot.ts(ss_start, ylab = 'count', main = 'SS, Start (blk), Final (bl), Mean (or)')
  lines(mcmc_output$ss[mcmc_vec_size, ], col = 'aquamarine', lwd = 2) #Final
  lines(colMeans(mcmc_output$ss), col = '#FF3300', lwd = 2) #mean

  #v. MEAN T1, T2
  plot.ts(colMeans(mcmc_output$ss[((mcmc_vec_size/2)+1):mcmc_vec_size, ]), ylab = 'count',
          main = 'SS - Mean T1 (prp), Mean T2 (or)', col = 'orange')
  lines(colMeans(mcmc_output$ss[1:(mcmc_vec_size/2), ]), col = 'blueviolet', lwd = 2)

  #SSE DATA
  } else {
    plot.ts(0, xlab = '', ylab = ''); plot.ts(0, xlab = '', ylab = '')
    plot.ts(0, xlab = '', ylab = ''); plot.ts(0, xlab = '', ylab = '')
  }
  
  #************************
  #ROW 2: MCMC TRACE PLOTS
  #************************
  
  #***************
  #r0
  plot.ts(r0_mcmc, ylab = "R0", #ylim=c(0, r0_lim),
          main = paste0("R0 Total. Simulated:", r0_sim),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = r0_sim, col = 'orange', lwd = 2)
  
  #a
  plot.ts(m1_mcmc, ylab = mod_par_names[1], #ylim=c(0, m1_lim),
          main = paste(mod_par_names[1], "trace.",
                       "Simulated: ", sim_vals$m1),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m1, col = 'red', lwd = 2)

  #***************
  #b
  plot.ts(m2_mcmc, ylab = mod_par_names[3], #ylim=c(0, max(m2_mcmc)),
          main = paste(mod_par_names[2], "trace.",
                       "Simulated: ", sim_vals$m2),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m2, col = 'blue', lwd = 2)

  #***************
  #c
  plot.ts(m3_mcmc,  ylab =  paste0(mod_par_names[3], ",sigma"), #ylim=c(0, m3_lim),
          main = paste(mod_par_names[3], "trace.", #m3_lim =  max(sim_vals$m3, max(m3_mcmc, na.rm = TRUE))
                       "Simulated: ", sim_vals$m3),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m3, col = 'green', lwd = 2)

  #***************
  #LOG LIKELIHOOD
  plot.ts(log_like_mcmc, ylab = "log likelihood",
          main = "Log Likelihood",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #**********************************************************
  #ROW 3:  HISTOGRAMS OF PARARMS (a, b, c, r0, loglike)
  #************************************************************
  
  #***********
  #HIST r0
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0', #ylab = 'Density',
       main = paste('R0, '," prior:", mr0_prior), #xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = r0_sim, col = 'orange', lwd = 2)
  #PRIOR
  if (FLAGS_LIST$PRIOR) {
    lines(xseq_r0, dr0e, type = 'l', lwd = 2) 
  }
  
  #HIST m1
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[1], #ylab = 'Density',
       main = paste(mod_par_names[1],
                    " prior:", m1_prior),
       #xlim=c(0, m1_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = sim_vals$m1, col = 'red', lwd = 2)

  #PRIOR PLOT
  if (FLAGS_LIST$PRIOR) {
    lines(xseq_a, dba, type = 'l', lwd = 2) 
  }

  #***********
  #HIST m2
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[2], #ylab = 'Density',
       main = paste(mod_par_names[2]),
                  #  " prior:", m2_prior),
      # xlim=c(0, m2_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = sim_vals$m2, col = 'blue', lwd = 2)

  #***********
  #Hist m3
  hist(m3_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[3], #ylab = 'Density',
       main = paste(mod_par_names[3],
                    " prior:", m3_prior),
      # xlim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = sim_vals$m3, col = 'green', lwd = 2)

  #PRIOR PLOT
  if (FLAGS_LIST$PRIOR) {
    lines(xseq_c, dgc, type = 'l', lwd = 2) 
  }

  #***********
  #HIST log_like_vec
  hist(log_like_mcmc, freq = FALSE, breaks = 100,
       xlab = 'Log likelihood', #ylab = 'Density',
       main = 'Log likelihood',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

  #************************************************
  #ROW 4: CUMULATIVE MEAN PLOTS
  #************************************************
  
  
  #r0 mean
  plot(seq_along(r0_mean), r0_mean, #ylim=c(0, r0_lim),
       xlab = 'Time', ylab = 'R0', main = paste0("R0 MCMC Mean, Simulated:", r0_sim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      # ylim = c(0, r0_lim),lwd = 1)
  abline(h = r0_sim, col = 'orange', lwd = 2)
  
  #m1 mean
  plot(seq_along(m1_mean), m1_mean,
       #ylim=c(0, m1_lim),
       xlab = 'Time', ylab =  mod_par_names[1],
       main = paste(mod_par_names[1], "MCMC mean, Start:", sim_vals$m1),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m1, col = 'red', lwd = 2)

  #m2 mean
  plot(seq_along(m2_mean), m2_mean,
       #ylim=c(0, m2_lim),
       xlab = 'Time', ylab = 'm2',
       main = paste(mod_par_names[2], "MCMC mean, Simulated:", sim_vals$m2),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  abline(h = sim_vals$m2, col = 'blue', lwd = 2)

  #m3 Mean
  plot(seq_along(m3_mean), m3_mean,
       xlab = 'Time', ylab = 'm3',
       main = paste(mod_par_names[3], "MCMC mean, Simulated:", sim_vals$m3),
       #ylim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  abline(h = sim_vals$m3, col = 'green', lwd = 2)

  #loglike mean
  plot(seq_along(cumsum(log_like_mcmc)/seq_along(log_like_mcmc)), cumsum(log_like_mcmc)/seq_along(log_like_mcmc),
       xlab = 'Time', ylab = 'log likelihood',
       main = "Log Likelihood Mean",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  
  #*****************
  #ROW 5: JOINT DISTRIBUTIONS/MARGINALS
  #********************

  #a vs r0
  plot(m1_mcmc, r0_mcmc,
       xlab = mod_par_names[1], ylab = 'R0',
       main = paste0(mod_par_names[1], ' vs R0'),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)

  #b*c vs r0
  plot(m2_mcmc*m3_mcmc, r0_mcmc,
       xlab = 'beta*gamma', ylab = 'R0',
       main = 'beta*gamma vs R0',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)

  #m1 vs m2
  plot(m1_mcmc, m2_mcmc,
       xlab = mod_par_names[1], ylab = mod_par_names[2],
       main = paste(mod_par_names[1], 'vs', mod_par_names[2]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)

  #m1 vs m3
  plot(m1_mcmc, m3_mcmc,
       xlab = mod_par_names[1], ylab = mod_par_names[3], main = paste(mod_par_names[1], 'vs', mod_par_names[3]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)

  #v. m2 vs m3
  plot(m2_mcmc, m3_mcmc,
       xlab = mod_par_names[2], ylab = mod_par_names[3], main = paste(mod_par_names[2], 'vs', mod_par_names[3]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       cex = 0.5)

  #********************
  #v. DATAFRAME       
  #********************
  df_results <- data.frame(
    #rep = mcmc_specs$seed_count,
    n_mcmc = n_mcmc,
    m1 = sim_vals$m1,
    m1_mc = round(mean(m1_mcmc), 2),
    m2 = sim_vals$m2,
    m2_mc = round(mean(m2_mcmc), 2),
    m3 = sim_vals$m3,
    m3_mc = round(mean(m3_mcmc), 2),
    R0 = r0_sim,
    R0_mc = round(mean(r0_mcmc), 2), 
    accept_rate_m1 = round(mcmc_output$list_accept_rates$accept_rate1, 2),
    a_rte_m2 = round(mcmc_output$list_accept_rates$accept_rate2, 2),
    a_rte_m3 = round(mcmc_output$list_accept_rates$accept_rate3, 2),
    a_es = effectiveSize(as.mcmc(m1_mcmc))[[1]],
    b_es = effectiveSize(as.mcmc(m2_mcmc))[[1]],
    c_es = effectiveSize(as.mcmc(m3_mcmc))[[1]],
    r0_es = effectiveSize(as.mcmc(r0_mcmc))[[1]])
    #time_elap = format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
    
    print(df_results)

  #return(df_results)

}
