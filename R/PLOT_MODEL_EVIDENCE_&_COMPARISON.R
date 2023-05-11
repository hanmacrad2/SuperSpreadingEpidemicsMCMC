#PLOTTING MODEL EVIDENCE RESULTS

#******************************************************************************
#* PLOTTING RESULTS FUNCTION
#*
#******************************************************************************
PLOT_BAYES_FACTORS <- function(bayes_factors, result_type = 'Bayes Factors via Harmonic Mean: Baseline vs SSEB Models. ',
                               data_type = 'Baseline', 
                               n_reps = 100, FLAG_RESULT_TYPE = list(log = FALSE)){
  
  #TITLE
  if (FLAG_RESULT_TYPE$log) {
    #axis_label = paste0(result_type, ' (log).')
    axis_label = paste0(result_type)
  } else  axis_label = paste0(result_type)
  
  #Title
  titleX = paste0(axis_label, data_type, ' data. ', n_reps, ' reps.')
  labelX =  'Bayes Factor'
  #PLOT
  #par(mfrow = c(2,1))
  boxplot(bayes_factors,
          ylab =labelX,
          main = axis_label)
  
  hist(bayes_factors, breaks = 100, freq = FALSE,
       xlab =labelX,
       main = axis_label)
  
}


#MODEL EVIDENCE RESULTS
PLOT_MODEL_EV_RESULTS <- function(posterior_results, model_type = 'SSEB',
                                  data_type = 'Baseline', 
                                  n_reps = 100, FLAG_RESULT_TYPE = list(phat = FALSE, post_prob = FALSE,
                                                                        log_model_ev = TRUE, log = FALSE)){
  
  #TITLE
  if(FLAG_RESULT_TYPE$phat) result_type = 'P hat, '
  if(FLAG_RESULT_TYPE$post_prob) result_type = 'Posterior model probability '
  if(FLAG_RESULT_TYPE$log_model_ev) result_type = 'Log Model Evidence '
  #LOG = TRUE
  if (FLAG_RESULT_TYPE$log) {
    axis_label = paste0(result_type, '(log), ', model_type, ' model. ')
    posterior_results = log(posterior_results)
  } else  axis_label = paste0(result_type, model_type, ' model. ')
  
  #Title
  titleX = paste0(axis_label, data_type, ' data. ', n_reps, ' reps.')
  
  #PLOT
  par(mfrow = c(2,1))
  boxplot(posterior_results,
          ylab = axis_label,
          main = titleX)
  
  hist(posterior_results, breaks = 50, freq = FALSE,
       xlab = axis_label,
       main = titleX)
  
}

#DATAFRAME OF POSTERIOR RESULTS
BOX_PLOT_POSTERIOR_PROBS <- function(list_vec_results = list(sseb = results1,
                                                          base = results2, ssnb = results3),
                                  data_type = data_type, model_ev_method = model_ev_method, titleX = '') { #Posterior Model Probabilities (Model evidence via Harmonic Mean). Data - Baseline Model
  
  #Set up
  title = paste0(titleX, 'Posterior Model Probabilities. ', data_type, model_ev_method) #' data. ', #' model evidence'
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  print(head(df_results))
  boxplot(df_results, main = title,
          col = c('red', 'green', 'blue'),
          cex.lab=1.3, cex.axis=1.3, cex.main=1.2, cex.sub=1.3)
  
  
}

#BOX PLOT RESULTS 
BOX_PLOT_MODEL_EV <- function(list_vec_results = list(importance_sampling = results1,
                                                      harmonic_mean = results2),
                                  data_type = data_type, model = model) { #Posterior Model Probabilities (Model evidence via Harmonic Mean). Data - Baseline Model
  
  #Set up
  title = paste0(model, ' Model Evidence. Importance sampling vs Harmonic Mean. ', data_type, ' data. ')
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  print(head(df_results))
  boxplot(df_results, main = title,
          col = c('red', 'green', 'blue'),
          cex.lab=1.3, cex.axis=1.3, cex.main=1.2, cex.sub=1.3)
  
  
}

#MODEL EVIDENCE RESULTS
BOX_PLOT_MODEL_EV_RESULTS <- function(list_mod_ev1, list_mod_ev2, list_mod_ev3, list_mod_ev4,
                                      model_type = 'SSEB', data_type = 'Baseline', 
                                      n_reps = 100, FLAG_RESULT_TYPE = list(phat = FALSE, post_prob = FALSE,
                                                                            log_model_ev = TRUE, log = FALSE)){
  
  #TITLE
  # if(FLAG_RESULT_TYPE$phat) result_type = 'P hat, '
  # if(FLAG_RESULT_TYPE$post_prob) result_type = 'Posterior model probability '
  # if(FLAG_RESULT_TYPE$log_model_ev) result_type = 'Log Model Evidence '
  # #LOG = TRUE
  # if (FLAG_RESULT_TYPE$log) {
  #   axis_label = paste0(result_type, '(log), ', model_type, ' model. ')
  #   posterior_results = log(posterior_results)
  # } else  axis_label = paste0(result_type, model_type, ' model. ')
  
  #Title
  #titleX = paste0(axis_label, data_type, ' data. ', n_reps, ' reps.')
  
  #PLOT
  par(mfrow = c(2,2)); ylabX = 'log( Model Evidence)'
  boxplot(list_mod_ev1,
          main = 'Harmonic Mean Model Evidence (log) Baseline model',
          ylab = ylabX)
  #2.
  boxplot(list_mod_ev2,
          main = 'Importance Sampling Model Evidence (log) Baseline model',
          ylab = ylabX)
  
  boxplot(list_mod_ev3,
          main = 'Harmonic Mean Model Evidence (log) SSEB model',
          ylab = ylabX)
  
  boxplot(list_mod_ev4,
          main = 'Importance Sampling Model Evidence (log) SSEB model',
          ylab = ylabX)
}