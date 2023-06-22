#INSPECT SD OF MODEL EVIDENCE ESTIMATE VS BETA

#FOLDER


SD_SSIB_MODEL_EV <- function(EPI_DATA, OUTER_FOLDER, n_repeats){
  
  #Beta 
  beta_list = c(0.01, 0.25, seq(0.05, 1, by = 0.05)) 
  #beta_list = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  sd_beta = vector(length = length(beta_list))
  
  for (i in 1:length(beta_list)){
    
    list_model_ev_ssib = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                      beta_ssib = beta_list[i],
                                 FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                     SSIB = TRUE, SSIR = FALSE))
    sd_beta[i] = sd(list_model_ev_ssib)
  }

  #PLOT
  plot(beta_list, sd_beta, main = 'Beta vs sd(Model Evidence of SSIB)',
       xlab = 'beta', ylab = 'sd(model evidence ssib)',
       pch = 16)
  
  return(list(beta_list = beta_list, sd_beta = sd_beta))
  
}

#PLOT
sd_output = SD_SSIB_MODEL_EV(data_ssib2, OUTER_FOLDER, n_repeats)

par(mfrow = c(1,1))
plot(sd_output$beta_list, sd_output$sd_beta, main = 'Beta vs sd(Model Evidence of SSIB)',
     xlab = 'beta', ylab = 'sd(model evidence ssib)',
     pch = 16)
