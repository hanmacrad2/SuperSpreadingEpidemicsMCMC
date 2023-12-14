#***************************
#PLOT MODEL EVIDENCE RESULTS

#*****************************
BAR_PLOT_POST_PROBS <- function(list_vec_results, title = '') {
  
  #PLOT MODEL PROBABILITIES
  
  #SETUP
  title = paste0('Posterior Model Probabilities. ', title)
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  
  bar_colors <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red') 
  
  barplot(t(df_results), beside = TRUE, col = bar_colors, ylim = c(0, 1),
          names.arg = colnames(df_results), main = title,
          xlab = 'Model', ylab = 'Posterior Probability',
          cex.lab=1.5, cex.axis=1.6, cex.main= 1.6, 
          cex.names = 1.6, cex.sub=1.6)
  legend('topright', legend = colnames(df_results), fill = bar_colors) #, cex = 1.1)
  
}


BOX_PLOT_POST_PROBS <- function(df_results, title = '') {
  
  # SETUP
  title <- paste0('Posterior Model Probabilities. ', title)
  
  box_colors <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red') 
  
  boxplot(df_results, col = box_colors, names = colnames(df_results),
          main = title, xlab = 'Model', ylab = 'Posterior Probability',
          cex.lab = 1.5, cex.axis = 1.6, cex.main = 1.6, 
          cex.names = 1.6, cex.sub = 1.6)
  
  legend('topright', legend = colnames(df_results), fill = box_colors)
}

BOX_PLOT_MOD_EV <- function(df_results, title = 'Model Evidence') {
  
  # SETUP
  #title <- paste0('Posterior Model Probabilities. ', title)
  
  box_colors <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red') 
  model_names = c('Baseline', 'SSE', 'SSI', 'SSEB', 'SSIB')
  
  boxplot(df_results, col = box_colors, names = model_names,
          main = title, xlab = 'Model', ylab = 'log Model Evidence',
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, 
          cex.names = 1.3, cex.sub = 1.3)
  
  #legend('topright', legend = colnames(df_results), fill = box_colors)
}

#PLOT TWO
BOX_PLOT_POSTERIOR_PROBS <- function(list_vec_results, titleX = '') { 
  
  'PLOT POSTERIOR MODEL PROBABILITY RESULTS'
  
  #Set up
  title = paste0(titleX, 'Posterior Model Probabilities. ', data_type, model_ev_method) #' data. ', #' model evidence'
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  print(head(df_results))
  boxplot(df_results, main = title,
          col = c('red', 'green', 'blue', 'orange', 'black'),
          ylim = c(0,1),
          cex.lab=1.3, cex.axis=1.3, cex.main=1.2, cex.sub=1.3)
  
  
}

#MODEL EVIDENCE
BOX_PLOT_MODEL_EV <- function(list_vec_results,
                              title = 'Model Evidence ',
                              data_type = '') { 
  
  'PLOT MODEL EVIDENCE RESULTS'
  #Set up
  title = paste0(title, data_type)
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  print(head(df_results))
  boxplot(df_results, main = title,
          col = c('red', 'green', 'blue', 'orange'),
          ylim = c((min(df_results)-0.5), (max(df_results)+0.5)),
          cex.lab=1.3, cex.axis=1.3, cex.main=1.2, cex.sub=1.3)
  
  
}

