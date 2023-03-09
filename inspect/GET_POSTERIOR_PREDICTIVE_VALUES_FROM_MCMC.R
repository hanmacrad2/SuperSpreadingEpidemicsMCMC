#GET POSTERIOR PREDICTIVE VALUES FROM MCMC

#PLOT POSTERIOR PREDICTIVE VALUES
PLOT_POSTERIOR_PRED_EPI_DATA <- function(true_epidemic_data, OUTER_FOLDER, true_R0 = 1.6, run = 1, n_repeats = 100,
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
  data_title = bquote("Posterior Predictive data. Base (orange), SSEB (green). " ~ bold(R[0] ~ .(true_R0)))
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = data_title, lwd = 3,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Parameters
  estimates_vec = c()
  
  if (FLAGS_LIST$BASE){
    #FOLDER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'BASE/run_', run) 
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      #READ MCMC SAMPLES
      mcmc_samples = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_base_', i ))
      mean_r0 = mean(mcmc_samples$r0_vec)
      posterior_pred_data = SIMULATE_BASELINE_EPIDEMIC(mean_r0)
      #PLOT
      lines(posterior_pred_data, col = 'orange')
      #SAVE
      saveRDS(posterior_pred_data, file = paste0(CURRENT_OUTPUT_FOLDER, '/post_pred_base_', i, '.rds' ))
    }
    #SAVE ESTIMATES
    
  } 
  
  if(FLAGS_LIST$SSEB){
    
    #FOLER 
    RESULTS_FOLDER = paste0(OUTER_FOLDER, 'SSEB/run_', run) 
    print('sseb')
    for (i in 1:n_repeats){
      
      print(paste0('i = ', i))
      mcmc_output = readRDS(file = paste0(RESULTS_FOLDER, '/mcmc_sseb_', i ))
      alpha_mean = mean(mcmc_output$alpha_vec);  beta_mean = mean(mcmc_output$beta_vec)
      gamma_mean = mean(mcmc_output$gamma_vec)
      #POSTERIOR PRED DATA
      posterior_pred_data = SIMULATE_EPI_SSEB(num_days = 50, alphaX = alpha_mean, betaX = beta_mean, gammaX = gamma_mean)
      #PLOT
      lines(posterior_pred_data, col = 'green')
      #SAVE
      saveRDS(posterior_pred_data, file = paste0(CURRENT_OUTPUT_FOLDER, '/post_pred_sseb_', i, '.rds' ))
      
    } 
  }
}

#Apply

    
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