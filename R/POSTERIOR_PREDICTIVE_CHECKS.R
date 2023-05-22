#GET POSTERIOR PREDICTIVE VALUES FROM MCMC
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"

#POSTERIOR PREDICTIVE PLOTS
POSTERIOR_PREDICTIVE_PLOTS <- function(matrix_sim_data, true_data){
  
  #SETUP
  par(mfrow = c(1, 2))
  num_days = length(true_data)
  
  #TEST STATS
  mean_est <- apply(matrix_sim_data, 2, mean) |> unlist()
  upper_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.975) |> unlist()
  lower_bounds <- apply(matrix_sim_data, 2, quantile, probs = 0.025) |> unlist()
  posterior_zig_zag = apply(m, 1, zigzag)
  
  #PLOTS 
  #HISTOGRAM OF ZIG-ZAG
  hist(posterior_zig_zag, xlim = c(0, max(max(posterior_zig_zag), posterior_zig_zag(true_data) * 2)),
       title = 'Zig-zag of data: sum(abs(diff(data)))')
  abline(v = quantile(true_data, probs = c(0.025, 0.975)), col = 'red')
  abline(v = zigzag(true_data))
  
  #PLOTS
  plot(1:num_days, upper_bounds, type = 'l', ylim = c(-5, 5))
  
  lines(1:num_days, rep(qnorm(p=0.025, sd = 2), num_days), type = 'l')
  
  lines(1:num_days, mean_est, type = 'l')
  
  lines(1:num_days, rep(qnorm(p=0.975, sd = 2), num_days), type = 'l')
  
  lines(1:num_days, lower_bounds, type = 'l')
  points(1:num_days, true_data, type = 'l', col = 'red')
}



#SAMPLE FROM BASELINE MODEL
#' @export
SAMPLE_BASELINE_MCMC <- function(mcmc_output, PLOT = FALSE){
  
  #SAMPLE
  r0_vec = mcmc_output$r0_vec
  #SAMPLE
  n_mcmc = length(r0_vec); sample_index = sample(1:n_mcmc, 1)
  #SAMPLE PARAMETERS
  r0 = r0_vec[sample_index]
  posterior_pred_data = SIMULATE_BASELINE_EPIDEMIC(r0)
  #PLOT
  if(PLOT)lines(posterior_pred_data, col = 'red')
  
  return(posterior_pred_data)
}

#SAMPLE FROM SSEB MODEL
#' @export
SAMPLE_SSEB_MCMC <- function(mcmc_output, PLOT = FALSE){
  
  #SAMPLE
  alpha_vec = mcmc_output$alpha_vec
  n_mcmc = length(alpha_vec)
  sample_index = sample(1:n_mcmc, 1)
  #SAMPLE PARAMETERS
  alpha = alpha_vec[sample_index]; beta = mcmc_output$beta_vec[sample_index]
  gamma = mcmc_output$beta_vec[sample_index]
  #POSTERIOR PRED DATA
  posterior_pred_data = SIMULATE_EPI_SSEB(num_days = 50, alphaX = alpha, betaX = beta, gammaX = gamma)
  
  #PLOT
  if(PLOT)lines(posterior_pred_data, col = 'green')
  
  return(posterior_pred_data)
}

#SAMPLE FROM SSEB MODEL
#' @export
SAMPLE_SSNB_MCMC <- function(mcmc_output, PLOT = FALSE){
  
  #SAMPLE
  n_mcmc = length(mcmc_output$ssnb_params_matrix[,1])
  sample_index = sample(1:n_mcmc, 1)
  #SAMPLE PARAMETERS
  kX = mcmc_output$ssnb_params_matrix[sample_index, 1]; R0X = mcmc_output$ssnb_params_matrix[sample_index, 2]
  #POSTERIOR PRED DATA
  posterior_pred_data = SIMULATE_EPI_SSNB(num_days = 50, R0 = R0X, k = kX)
  
  #PLOT
  if(PLOT)lines(posterior_pred_data, col = 'blue')
  
  return(posterior_pred_data)
}

#PLOT POSTERIOR PREDICTIVE DATA
#' @export
PLOT_POSTERIOR_PRED_EPI_DATA <- function(true_epidemic_data, OUTER_FOLDER, true_R0 = 1.6,
                                         run = 1, n_repeats = 100, n_sample_repeats = 5, burn_in = 2000,
                                         FLAGS_MODELS = list(BASE = TRUE, SSEB = TRUE,
                                                             SSNB = TRUE, SSIB = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #FOLDER
  #CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run) 
  #create_folder(CURRENT_OUTPUT_FOLDER)
  
  #PLOT TRUE
  par(mfrow = c(1,1))
  count = 0
  data_title = bquote("Posterior Predictive data. Baseline model (red), SSEB model (green), SSNB model (blue). " ~ bold(R[0]: ~ .(true_R0)))
  plot.ts(true_epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          ylim = c(0,100),
          main = data_title, lwd = 2,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #1.SSEB PREDICTIVE DATA
  if(FLAGS_MODELS$SSEB){
    
    #FOLER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'SSEB/run_', run) 
    print('sseb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #MCMC OUTPUT
      mcmc_output = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_sseb_', i ))
      
      for (j in 1:n_sample_repeats){
        SAMPLE_SSEB_MCMC(mcmc_output, PLOT = TRUE)
        count = count+ 1
      }
    } 
  }
  
  if(FLAGS_MODELS$SSNB){
    
    #FOLER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'SSNB/run_', run) 
    print('ssnb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #MCMC OUTPUT
      print(paste0(RESULTS_FOLDER, '/mcmc_', i ))
      mcmc_output = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_', i ,'.rds'))
      
      for (j in 1:n_sample_repeats){
        SAMPLE_SSNB_MCMC(mcmc_output, PLOT = TRUE)
        count = count+ 1
      }
    } 
  }
  
  #2. BASE PREDICTIVE DATA
  if (FLAGS_MODELS$BASE){
    #FOLDER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'BASE/run_', run) 
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #READ MCMC SAMPLES
      mcmc_samples = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_base_', i ))
      r0_vec = mcmc_output$r0_vec
      
      for (j in 1:n_sample_repeats){
        SAMPLE_BASELINE_MCMC(mcmc_output, PLOT = TRUE)
        count = count + 1
      }
    }
  } 
  #PLOT ON TOP
  lines(true_epidemic_data, #xlab = 'Time', ylab = 'Daily Infections count', main = data_title
        lwd = 2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}

#Apply
#PLOT_POSTERIOR_PRED_EPI_DATA(data_baseI, OUTER_FOLDER)

#POSTERIOR PREDICTIVE CHECKING
PLOT_POSTERIOR_PRED_EPI_DATA <- function(true_epidemic_data, OUTER_FOLDER,
                                         run = 1, num_reps = 100,
                                         MODELS = list(BASELINE = 'BASELINE', SSEB = 'SSEB',
                                                       SSNB = 'SSNB', SSIB = 'SSIB', SSIR = 'SSIR'))
                                         FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE,
                                                             SSNB = FALSE, SSIB = FALSE, SSIR = FALSE)){
  #STORE MATRIX DATA
  matrix_data <- matrix(0, nrow = num_reps, ncol = num_days)
  
   if (FLAGS_MODELS$BASELINE){
    #FOLDER
    model_type = 'BASELINE'
    RESULTS_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run) 
    
    for (i in 1:n_repeats){
      
      #READ MCMC SAMPLES
      mcmc_output = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_base_', i ))
      matrix_data[i] = SAMPLE_BASELINE_MCMC(mcmc_output)

    }
   } 
  
  #POSTERIOR PREDICTIVE PLOTS
  

}