#***************************
#PLOT MODEL EVIDENCE RESULTS
#*****************************

#************************************
# SINGLE MODEL COMPARISON RUNS
#************************************
BAR_PLOT_POST_PROBS <- function(list_vec_results, title = '') {
  
  #PLOT MODEL PROBABILITIES
  
  #SETUP
  title = paste0(title, ' Posterior Model Probabilities. ')
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
  title <- paste0(title, ' Posterior Model Probabilities. ')
  
  box_colors <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red') 
  
  boxplot(df_results, col = box_colors, names = colnames(df_results),
          main = title, xlab = 'Model', ylab = 'Posterior Probability',
          cex.lab = 1.5, cex.axis = 1.6, cex.main = 1.6, 
          cex.names = 1.6, cex.sub = 1.6)
  
  legend('topright', legend = colnames(df_results), fill = box_colors)
}

#************************************
# MULTIPLE MODEL COMPARISON RUNS
#************************************
BOX_PLOT_MOD_EV_MULT <- function(df_results, title = 'Model Evidence') {
  
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

#POSTERIOR PROBS
BOX_PLOT_POST_PROBS_MULT <- function(df_results, title = '') {
  
  # SETUP
  title <- paste0(title, ' - Posterior Model Probabilities')
  
  box_colors <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red') 
  model_names = c('Baseline', 'SSE', 'SSI', 'SSEB', 'SSIB')
  
  boxplot(df_results, col = box_colors, names = model_names,
          main = title, xlab = 'Model', ylab = 'Posterior Model Probability',
          cex.lab = 1.5, cex.axis = 1.6, cex.main = 1.6, 
          cex.names = 1.6, cex.sub = 1.6)
  #legend('topright', legend = colnames(df_results), fill = box_colors)
}

