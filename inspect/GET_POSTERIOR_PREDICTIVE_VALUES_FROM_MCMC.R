#GET POSTERIOR PREDICTIVE VALUES FROM MCMC

#PLOT POSTERIOR PREDICTIVE VALUES
PLOT_POSTERIOR_PRED_EPI_DATA <- function(true_epidemic_data, OUTER_FOLDER, true_R0 = 1.6, run = 1, n_repeats = 100, burn_in = 2000,
                                FLAGS_MODELS = list(BASE = TRUE, SSEB = TRUE,
                                                    SSIB = FALSE, SSIC = FALSE)){
  'For a given epidemic dataset and model. 
  Get importance sampling estimate of model evidence. 
  1. Run mcmc 2. Get estimate'
  
  #FOLDER
  #CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run) 
  #create_folder(CURRENT_OUTPUT_FOLDER)
  
  #PLOT TRUE
  par(mfrow = c(1,1))
  data_title = bquote("Posterior Predictive data. Baseline model (orange), SSEB model (green). " ~ bold(R[0]: ~ .(true_R0)))
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
      alpha_vec = mcmc_output$alpha_vec
      #SAMPLE
      n_mcmc = length(alpha_vec); sample_index = sample(burn_in:n_mcmc, 1)
      #SAMPLE PARAMETERS
      alpha = alpha_vec[sample_index]; beta = mcmc_output$beta_vec[sample_index]
      gamma = mcmc_output$beta_vec[sample_index]
      #POSTERIOR PRED DATA
      posterior_pred_data = SIMULATE_EPI_SSEB(num_days = 50, alphaX = alpha, betaX = beta, gammaX = gamma)
      #PLOT
      lines(posterior_pred_data, col = 'green')
      #SAVE
      saveRDS(posterior_pred_data, file = paste0(RESULTS_FOLDER, '/post_pred_sseb_', i, '.rds' ))
      
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
      #SAMPLE
      n_mcmc = length(r0_vec); sample_index = sample(burn_in:n_mcmc, 1)
      #SAMPLE PARAMETERS
      r0 = r0_vec[sample_index]; 
      posterior_pred_data = SIMULATE_BASELINE_EPIDEMIC(r0)
      #PLOT
      lines(posterior_pred_data, col = 'orange')
      #SAVE
      saveRDS(posterior_pred_data, file = paste0(RESULTS_FOLDER, '/post_pred_base_', i, '.rds' ))
    }
    #SAVE ESTIMATES
    
  } 
  #PLOT ON TOP
  lines(true_epidemic_data, #xlab = 'Time', ylab = 'Daily Infections count', main = data_title
        lwd = 2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}

#Apply
PLOT_POSTERIOR_PRED_EPI_DATA
    
  # } else if (FLAGS_LIST$SSIB){
  #   
  #   for (i in 1:n_repeats){
  #     
  #     print(paste0('i = ', i))
  #     mcmc_output = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_ssib_', i ))
  #     mcmc_samples =  matrix(c(mcmc_output$a_vec, mcmc_output$b_vec, mcmc_output$c_vec), ncol = 3)
  #     
  #     #GET PHAT ESTIMATE OF MODEL EVIDENCE
  #     phat_estimate = GET_LOG_P_HAT(mcmc_samples, epidemic_data) 
  #     estimates_vec[i] = phat_estimate
  #     print(estimates_vec)
  #     
  #   }
  #   
  #   #SAVE ESTIMATES
  #   saveRDS(estimates_vec, file = paste0(CURRENT_OUTPUT_FOLDER, '/phat_ests_ssib_', run, '.rds' ))
  #   
  # }