#PLOT SSIR GRID

#PLOT ETA ACROS DAYS
PLOT_ETA <- function(epidemic_data, eta_matrix, true_eta_vec,
                     seedX, eta_start = 1,
                     eta_reps = 16){
  
  #PLOT
  plot.new()
  par(mfrow=c(4,5))
  colours <- rainbow(eta_reps + 1); count = 1
  
  #PLOT INFECTIONS
  inf_title = paste0(" Data. Seed = ", seedX) 
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = inf_title,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #PLOT ETA
  plot.ts(true_eta_vec, xlab = 'Time', ylab = 'Daily Eta',
          main = "Eta: Epidemic Simulation",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #PLOT CREDIBLE INTERVAL
  ETA_CREDIBLE_INTERVALS(eta_matrix, true_eta_vec)
  
  for(day in seq(eta_start, eta_start + eta_reps, by = 1)){
    
    print(day)
    eta_dfX = eta_matrix[, day]
    
    plot.ts(eta_dfX,  ylab = paste0('Eta day ', day), #ylim = m3_lim,
            main = paste0('Eta day ', day),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h = true_eta_vec[day], col = colours[count], lwd = 2)
  count = count + 1
  
  }
}

#*****************************
#PLOT ETA CREDIBLE INTERVALS
ETA_CREDIBLE_INTERVALS <- function(eta_matrix, eta_true, pchX = 16,
                                   lwdX = 1){
  
  #Create a vector of means across columns
  eta_means = colMeans(eta_matrix)
  #Upper & lower limits
  ci = get_ci_matrix(eta_matrix) 
  
  #Plot
  plotCI(seq_along(eta_true), eta_means, ui = ci$vec_upper, li = ci$vec_lower,
         xlab = 'Day of Epidemic', ylab = 'Eta',
         main = 'Eta MCMC across the Epidemic 
       Posterior Mean & 95 % CI. Red (True/Simulated) ', lwd = lwdX, pch = 16) #xlim = c(min(vec_alpha), max(vec_alpha)))
  points(eta_true, col = 'red', lwd = lwdX, pch = pchX)
  lines(eta_true, col = 'red', lwd = lwdX)
  
  
}

#' Grid Plot of MCMC Individual R0; nu model
#'
#'Grid Plot of MCMC Results for  Individual R0; nu model
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
#' @param priors_list A list of prior parameters used
#' \itemize{
#'   \item \code{"a_prior_exp"} - rate of exponential prior on a, default rate of 1
#'   \item \code{"b_prior_exp"} - rate of exponential prior on b, default rate of 1
#'   \item \code{"b_prior_ga"} -  shape and scale of gamma distriubtion prior on b, defaults Ga(10, 0.02)
#'   \item \code{"c_prior_exp"}  - rate of exponential prior on c, default rate of 0.1
#'   \item \code{"c_prior_ga"} -  shape and scale of gamma distriubtion prior on c, defaults Ga(10, 1)
#' }
#' @param FLAGS_LIST A list of boolean variables for switching on/off certain functionality
#' \itemize{
#' \item \code{"BURN_IN"}  - Burn-in applied to mcmc samples if TRUE of size \code{"mcmc_specs$burn_in_pc"}
#'   \item \code{"THIN"}  - Return a thinned mcmc sample if TRUE, reduced by a factor of \code{"thinning_factor"}
#'   \item \code{"DATA_AUG"} - Data Augmentation was implemented as part of the SSI model if TRUE
#'   \item \code{"ADAPTIVE"} - Adaptive Algorithm applied to MCMC samples if TRUE
#'   \item \code{"MULTI_ALG"}  -  Multi Adaptive Shaping Algorithm applied to MCMC samples if TRUE
#'   \item \code{"PRIOR"}  - Plot prior distributions in grid if TRUE
#'   \item \code{"B_PRIOR_GAMMA"}  - A Gamma prior on b if TRUE, otherwise exponential
#'   \item \code{"C_PRIOR_GAMMA"}  - A Gamma prior on c if TRUE, otherwise exponential
#' }
#' @return Dataframe of results including mcmc sample starting points, mcmc sample final means, acceptance rates, mcmc effective sizes and the mcmc sampler run time
#' @export
#'
#'@author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' mcmc_plot_inputs = list(n_mcmc = 500000,  burn_in_pc = 0.05, mod_start_points = list(m1 = 0.72, m2 = 0.0038, m3 = 22),
#' mod_par_names = c('a', 'b', 'c'), model_type = 'SSI', seed_count = seed_count, thinning_factor = 10)
#'
#' df_mcmc_results = PLOT_SS_MCMC_GRID(epidemic_data, mcmc_output) 

#PLOT MCMC GRID
PLOT_SSI_MCMC_GRID <- function(epidemic_data, mcmc_output, eta_sim, seed_count,
                                log_like_sim, n_mcmc,
                                sim_vals = list(m1 = 1.6, m2 = 0.16),
                               mod_par_names = c('R0', 'k', 'eta'),
                               mcmc_specs = list(model_type = 'SSI',
                                                 mod_start_points = list(m1 = 1.2, m2 = 0.16),
                                                 burn_in_pc = 0.2,
                                                 eta_time_point = 1), #80 28
                               priors_list = list(k_prior = c(1,0), alpha_prior = c(1,0)),
                               FLAGS_LIST = list(PRIOR = TRUE,ADAPTIVE = FALSE, MULTI_ALG = TRUE,
                                                 MARGINALS = FALSE)){
  
  #FIX
  print(sim_vals)
  eta_sim_val = 0
  
  #PLOT
  plot.new()
  par(mfrow=c(4,4))
  
  #PRIORS mod_par_names = c('R0', 'k'),
  PRIORS_USED =  GET_PRIORS_USED()
  list_priors = GET_LIST_PRIORS_SSI() 
  
  #PRIOR R0
  if(PRIORS_USED$SSI$R0$EXP){
    mr0_prior = paste0('exp(', list_priors$r0[1], ')')
    xseq_r0 = seq(0, 3, length.out = 500)
    dr0e = dexp(xseq_r0, list_priors$r0[1])
  }
  
  #PRIOR k 
  if(PRIORS_USED$SSI$k$EXP){
    m2_prior = paste0('exp(', list_priors$k[1], ')')
    x2 = seq(0, 3, length.out = 500)
    d2 = dexp(x2, list_priors$k[1])
  }
  
  #EXTRACT MCMC SAMPLES
  log_like_mcmc = mcmc_output$log_like_vec; log_like_mcmc = unlist(log_like_mcmc)
  
  if (FLAGS_LIST$MULTI_ALG){
    m1_mcmc = mcmc_output$ssir_params_matrix[,1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output$ssir_params_matrix[,2]; m2_mcmc = unlist(m2_mcmc); m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
    m3_mcmc = mcmc_output$eta_matrix[, mcmc_specs$eta_time_point]; m3_mcmc = unlist(m3_mcmc); m3_mcmc = m3_mcmc[!is.na(m3_mcmc)]
    sigma_eta_X = mcmc_output$sigma_eta_matrix[, mcmc_specs$eta_time_point]; sigma_eta_X = unlist(sigma_eta_X); sigma_eta_X = sigma_eta_X[!is.na(sigma_eta_X)]
    #eta_sim_val = eta_sim[mcmc_specs$eta_time_point]
    #m3_mcmc = mcmc_output$x_matrix[,3]; m3_mcmc = unlist(m3_mcmc); m3_mcmc = m3_mcmc[!is.na(m3_mcmc)]
    
  } else {
    m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc); m1_mcmc = m1_mcmc[!is.na(m1_mcmc)]
    m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc);  m2_mcmc = m2_mcmc[!is.na(m2_mcmc)]
    m3_mcmc = mcmc_output[3]; m3_mcmc = unlist(m3_mcmc);  m3_mcmc = m3_mcmc[!is.na(m3_mcmc)]
  }
  
  
  #LIMITS
  m1_min =  min(mcmc_specs$mod_start_points$m1, min(m1_mcmc, na.rm = TRUE));  m1_max =  max(mcmc_specs$mod_start_points$m1, max(m1_mcmc, na.rm = TRUE))
  m2_min = min(mcmc_specs$mod_start_points$m2, min(m2_mcmc, na.rm = TRUE)); m2_max = max(mcmc_specs$mod_start_points$m2, max(m2_mcmc, na.rm = TRUE))
  m3_min = min(eta_sim_val, min(m3_mcmc, na.rm = TRUE)); m3_max = max(eta_sim_val, max(m3_mcmc, na.rm = TRUE))
  minll = min(min(log_like_mcmc, na.rm = TRUE), 0) #log_like_sim); 
  maxll = max(max(log_like_mcmc,  na.rm = TRUE), 0 ) #log_like_sim)
  #LIMITS
  m1_lim = c(m1_min, m1_max);  m2_lim = c(m2_min, m2_max);  m3_lim = c(m3_min, m3_max); m4_lim = c(minll, maxll)
  title_eta = paste(mod_par_names[3], "MCMC.", " Day: ", mcmc_specs$eta_time_point)
  #Priors
  m2_prior =  paste0('exp(', priors_list$k_prior[1], ')')
  
  
  #******************************************************************
  #* PLOTS *
  #******************************************************************
  
  #*****************
  #ROW 1: DATA INFECTIONS + JOINT DISTRIBUTIONS/MARGINALS
  #********************
  
  #i. TOTAL INFECTIONS
  inf_tite = paste0(mcmc_specs$model_type, " Data") #. Seed = ", seed_count) 
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = inf_tite,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  if(FLAGS_LIST$MARGINALS){
    
    #ALPHA VS K
    plot(m1_mcmc, m2_mcmc,
         xlab = mod_par_names[1], ylab = mod_par_names[2],
         main = paste0(mod_par_names[1], ' vs ', mod_par_names[2]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex = 0.5)
    
    #ETA vs K
    plot(m2_mcmc, m3_mcmc,
         xlab = mod_par_names[2], ylab = mod_par_names[3],
         main = paste0(mod_par_names[2], ' vs ', mod_par_names[3],
                       ' Day ', mcmc_specs$eta_time_point),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex = 0.5)
    
    #ETA VS ALPHA
    plot(m1_mcmc, m3_mcmc,
         xlab = mod_par_names[1], ylab = mod_par_names[3],
         main = paste0(mod_par_names[1], ' vs ', mod_par_names[3],
                       ' Day ', mcmc_specs$eta_time_point),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex = 0.5)
    
    
  } else {
    plot.ts(0, xlab = '', ylab = ''); plot.ts(0, xlab = '', ylab = '')
    plot.ts(0, xlab = '', ylab = '')
  }

  #************************
  #ROW 2: MCMC TRACE PLOTS
  #************************

  #R0
  r0_title = bquote(bold(R[0] ~ "MCMC, sim: " ~ .(sim_vals$m1)))
  #alpha_title = bquote(bold(alpha ~ "; Mean of" ~ eta ~ "~ Ga(). Simulated: " ~ .(sim_vals$m1)))
  
  if (!FLAGS_LIST$ADAPTIVE){
    plot.ts(m1_mcmc, ylab = mod_par_names[1], ylim= m1_lim, #bquote("Hello" ~ r[xy] == .(cor) ~ "and" ~ B^2)
            main = paste(mod_par_names[1], "MCMC",
                         "Simulated: ", sim_vals$m1),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h = sim_vals$m1, col = 'red', lwd = 2) #True = green
  
    } else {
    sig1 = mcmc_output$sigma$sigma1_vec
    plot.ts(m1_mcmc, ylab = paste0(mod_par_names[1], ",sigma"), #ylim=c(min(min(sig1),min(m1_mcmc)), max(m1_mcmc)),
            main = paste(mod_par_names[1], "MCMC",
                         "Start: ", sim_vals$m1, ', Sigma (red)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma1_vec, col = 'red')
  }
  
  #***************
  #k
  if (!FLAGS_LIST$ADAPTIVE){
    plot.ts(m2_mcmc, ylab = mod_par_names[3], ylim= m2_lim,
            main = paste(mod_par_names[2], "MCMC",
                         "Simulated: ", sim_vals$m2),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h = sim_vals$m2, col = 'blue', lwd = 2) #True = green
  } else {
    plot.ts(m2_mcmc, ylab = paste0(mod_par_names[2], ",sigma"), #ylim=c(0, m2_lim),
            main = paste(mod_par_names[2], "MCMC",
                         "Start: ", sim_vals$m2, ', Sigma (blue)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma2_vec, col = 'blue')
  }
  
  #***************
  #ETA
  if (!FLAGS_LIST$ADAPTIVE){
    plot.ts(m3_mcmc,  ylab = mod_par_names[3], ylim = m3_lim,
            main = title_eta,
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #abline(h = eta_sim_val, col = 'green', lwd = 2)
  } else {
    plot.ts(m3_mcmc,  ylab =  paste0(mod_par_names[3], ",sigma"), #ylim=c(0, m3_lim),
            main = paste(mod_par_names[3], "MCMC",
                         "Start: ", mcmc_specs$mod_start_points$m3, ', Sigma (green)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    lines(mcmc_output$sigma$sigma3_vec, col = 'green')
  }
  
  #***************
  #LOG LIKELIHOOD
  plot.ts(log_like_mcmc, ylab = "log likelihood",# ylim= m4_lim,
          main = paste("Log Likelihood. N MCMC:", n_mcmc, ". Burn-in:", mcmc_specs$burn_in_pc),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #abline(h = log_like_sim, col = 'orange', lwd = 2)
  
  #**********************************************************
  #ROW 2:  HISTOGRAMS OF PARARMS (a, b, c, r0, loglike)
  #************************************************************
  
  #***********
  #HIST m1
  #alpha_titleII = bquote(bold(alpha ~ "Prior: exp("~.(priors_list$alpha_prior[1])~")"))
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[1], #ylab = 'Density',
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
  #ktitle = bquote(bold(k ~ "Prior: exp("~.(priors_list$k_prior[1])~")"))
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[2], #ylab = 'Density',
       main = paste(mod_par_names[2],
                    " prior:", m2_prior),
       xlim= m2_lim,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v =sim_vals$m2, col = 'blue', lwd = 2)

  #PRIOR PLOT 
  if(FLAGS_LIST$PRIOR){
    lines(x2, d2, type = 'l', lwd = 2) 
  }
  
  
  #***********
  #Hist M3 (ETA)
  hist(m3_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[3], #ylab = 'Density',
       main = title_eta, xlim = m3_lim,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #                   " prior:", m3_prior),
  #abline(v = eta_sim_val, col = 'green', lwd = 2)
  
  #***********
  #HIST LOG_LIKE_VEC
  hist(log_like_mcmc, freq = FALSE, breaks = 100, #xlim= m4_lim,
       xlab = 'Log likelihood', #ylab = 'Density',
       main = 'Log likelihood',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #abline(v = log_like_sim, col = 'orange', lwd = 2)
  
  #************************************************
  #ROW 3: CUMULATIVE MEAN PLOTS
  #************************************************
  
  #m1 mean
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc)
  plot(seq_along(m1_mean), m1_mean,
       ylim= m1_lim,
       xlab = 'Time', ylab =  mod_par_names[1],
       main = paste(mod_par_names[1], "MCMC mean, Start:", mcmc_specs$mod_start_points$m1),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = sim_vals$m1, col = 'red', lwd = 2)
  
  #m2 mean
  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  plot(seq_along(m2_mean), m2_mean,
       ylim= m2_lim,
       xlab = 'Time', ylab = mod_par_names[2],
       main = paste(mod_par_names[2], "MCMC mean, Start:", mcmc_specs$mod_start_points$m2),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       lwd = 1)
  abline(h = sim_vals$m2, col = 'blue', lwd = 2)
  
  #ETA MEAN + CREDIBLE INTERVALS
  #!!!!!!!!!!!!!
  #INCLUDE BACK IN!!!!!!!!!!!!
  #ETA_CREDIBLE_INTERVALS(mcmc_output$eta_matrix, eta_sim, lwdX = 1)
  
  # #SIGMA ETA
  # plot.ts(sigma_eta_X,  ylab = paste0('sigma ', mod_par_names[3]), 
  #         main = paste0('Sigma ', title_eta),
  #         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #********************
  #v. DATAFRAME: RESULTS
  #********************
  
  df_results <- data.frame(
    #rep = seed_count,
    n_mcmc = n_mcmc,
    #mcmc_vec_size = mcmc_vec_size,
    alpha_sim = sim_vals$m1,
    #alpha_start = mcmc_specs$mod_start_points$m1[[1]],
    alpha_mean_mcmc = round(mean(m1_mcmc), 2), #round(mean(m1_mcmc[(mcmc_vec_size/2): mcmc_vec_size]), 2),
    k_sim = sim_vals$m2, 
    #k_start =  mcmc_specs$mod_start_points$m2[[1]],
    k_mean_mcmc = round(mean(m2_mcmc), 2),
    #eta_sim = round(eta_sim_val, 2),
    #eta_mean_mcmc = round(mean(m3_mcmc, na.rm = TRUE), 2),
    #log_like_sim = round(log_like_sim, 2),
    accept_rate = round(mcmc_output$accept_rate, 2),
    a_rte_d_aug = round(mcmc_output$accept_rate_da, 2),
    #alpha_es = round(effectiveSize(as.mcmc(m1_mcmc))[[1]], 2),
    #k_es = round(effectiveSize(as.mcmc(m2_mcmc))[[1]], 2),
    #eta_es = round(effectiveSize(as.mcmc(m3_mcmc))[[1]], 2),
    time_elap = mcmc_output$time_elap) #format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
  
  print(df_results)
  
  return(df_results)
  
}
