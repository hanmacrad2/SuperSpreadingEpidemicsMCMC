#GET POSTERIOR PREDICTIVE VALUES FROM MCMC
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"

#SAMPLE FROM BASELINE MODEL
SAMPLE_BASELINE_MCMC <- function(mcmc_output){
  
  #SAMPLE
  r0_vec = mcmc_output$r0_vec
  #SAMPLE
  n_mcmc = length(r0_vec); sample_index = sample(1:n_mcmc, 1)
  #SAMPLE PARAMETERS
  r0 = r0_vec[sample_index]
  posterior_pred_data = SIMULATE_BASELINE_EPIDEMIC(r0)
  #PLOT
  #lines(posterior_pred_data, col = 'red')
  
  return(posterior_pred_data)
}

#SAMPLE FROM SSEB MODEL
SAMPLE_SSEB_MCMC <- function(mcmc_output){
  
  #SAMPLE
  n_mcmc = length(mcmc_output$alpha_vec)
  sample_index = sample(1:n_mcmc, 1)
  #SAMPLE PARAMETERS
  alpha = mcmc_output$alpha_vec[sample_index]; beta = mcmc_output$beta_vec[sample_index]
  gamma = mcmc_output$beta_vec[sample_index]
  #POSTERIOR PRED DATA
  posterior_pred_data = SIMULATE_EPI_SSEB(alphaX = alpha, betaX = beta, gammaX = gamma)
  #PLOT
  #lines(posterior_pred_data, col = 'green')
  
  return(posterior_pred_data)
}

#SAMPLE FROM SSEB MODEL
SAMPLE_SSNB_MCMC <- function(mcmc_output){
  
  #SAMPLE
  n_mcmc = length(mcmc_output$ssnb_params_matrix[,1])
  sample_index = sample(1:n_mcmc, 1)
  #SAMPLE PARAMETERS
  kX = mcmc_output$ssnb_params_matrix[sample_index, 1]; R0X = mcmc_output$ssnb_params_matrix[sample_index, 2]
  #POSTERIOR PRED DATA
  posterior_pred_data = SIMULATE_EPI_SSNB(num_days = 50, R0 = R0X, k = kX)
  #PLOT
  #lines(posterior_pred_data, col = 'blue')
  
  return(posterior_pred_data)
}

# PLOT PRED DATA
PLOT_CI_PRED_DATA <- function(matrix_data_pred, colorX){
  
  #Matrix of data (confidence intervals)
  ci_matrix = get_ci_matrix(matrix_data_pred)
  #Plot
  lines(ci_matrix$vec_upper, type = "l", lty = "dashed", col = colorX, lwd = 2)
  lines(ci_matrix$vec_lower, type = "l", lty = "dashed", col =colorX, lwd = 2)
}


#PLOT POSTERIOR PREDICTIVE DATA
PLOT_POSTERIOR_PRED_EPI_DATA <- function(true_epidemic_data, OUTER_FOLDER, run, true_R0 = 1.6, num_epi_days = 50, 
                                         n_repeats = 100, n_sample_repeats = 5, upper_lim = 20, 
                                FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE,
                                                    SSNB = FALSE, SSIB = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #FOLDER
  
  #PLOT TRUE
  #par(mfrow = c(2,1))
  matrix_data_pred = matrix(0, n_repeats*n_sample_repeats, num_epi_days)
  count = 0
  #data_title = bquote("Posterior Predictive data (95% Credible Intervals), SSEB (green; true model), Baseline(red), SSNB (blue). " ~ bold(R[0]: ~ .(true_R0)))
  data_title = bquote("Posterior Predictive data (95% Credible Intervals), SSEB (green; true model). " ~ bold(R[0]: ~ .(true_R0)))
  plot.ts(true_epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          ylim = c(0, upper_lim),
          main = data_title, lwd = 2,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.1, cex.sub=1.5)
  
  #1.SSEB PREDICTIVE DATA
  if(FLAGS_MODELS$SSEB){
    
    #FOLER
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'SSEB/run_', run) 
    print('sseb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #MCMC OUTPUT
      mcmc_output = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_sseb_', i, '.rds'))
      
      for (j in 1:n_sample_repeats){
        data_pred = SAMPLE_SSEB_MCMC(mcmc_output)
        matrix_data_pred[(i*j),] = data_pred
        #count = count+ 1
      }
    } 
    #Matrix of data (confidence intervals)
    #print(matrix_data_pred)
    PLOT_CI_PRED_DATA(matrix_data_pred, 'green')
  }
  
  if(FLAGS_MODELS$SSNB){
    
    #FOLER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'SSNB/run_', run) 
    print('ssnb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #MCMC OUTPUT
      #print(paste0(RESULTS_FOLDER, '/mcmc_', i ))
      mcmc_output = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_ssnb_', i,'.rds'))
      
      for (j in 1:n_sample_repeats){
        data_pred = SAMPLE_SSNB_MCMC(mcmc_output)
        matrix_data_pred[(i*j),] = data_pred
        #count = count+ 1
      }
    }
    #Matrix of data (confidence intervals)
    PLOT_CI_PRED_DATA(matrix_data_pred, 'blue')
  }
  
  #2. BASE PREDICTIVE DATA
  if (FLAGS_MODELS$BASE){
    #FOLDER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'BASELINE/run_', run) 
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #READ MCMC SAMPLES
      mcmc_samples = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_baseline_', i,'.rds'))
      r0_vec = mcmc_output$r0_vec
      
      for (j in 1:n_sample_repeats){
        data_pred = SAMPLE_BASELINE_MCMC(mcmc_output)
        matrix_data_pred[(i*j),] = data_pred
        #count = count + 1
      }
    }
    #Matrix of data (confidence intervals)
    PLOT_CI_PRED_DATA(matrix_data_pred, 'red')
  } 
  
  #PLOT ON TOP
  lines(true_epidemic_data, #xlab = 'Time', ylab = 'Daily Infections count', main = data_title
        lwd = 2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}

#Apply
run = 2
PLOT_POSTERIOR_PRED_EPI_DATA(data_sseb, OUTER_FOLDER, run)
