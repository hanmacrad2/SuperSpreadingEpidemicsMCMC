#***************************
#PLOT MODEL EVIDENCE RESULTS
#*****************************
#library(vioplot)

#************************************
# SINGLE MODEL COMPARISON RUNS
#************************************
VIOLIN_PLOT_MODEL_COMPARISON <- function(df_results, RESULTS_FOLDER,
                                         cex = 2.0,
                                         plot_width = 12.5, plot_height = 10.0,
                                         MODEL_NAMES = c('Baseline', 'SSE', 'SSI', 'SSEB', 'SSIB'),
                                         MODEL_COLORS = c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', '#DC143C'),
                                         title = '', PLOT_VIOLIN = TRUE, PDF = TRUE) {
  
  # SETUP
  #old_par <- par(cex = cex)
  main_title = bquote(paste(.(model), " simulated data - Posterior Probabilities of the models"))
  xlab = bquote('Model')
  ylab = bquote('Posterior Model Probability')
  #title <- paste0(title, ' - Posterior Model Probabilities')
  
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, '/', model, '/plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0(model, '_MODEL_COMPARISON_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
  }
  
  #VIOLIN PLOT
  if(PLOT_VIOLIN){
    names(df_results) <- MODEL_NAMES
    par(cex = cex)
    vioplot(df_results, 
            col = MODEL_COLORS, names = MODEL_NAMES,
            main = main_title,
            xlab = xlab, 
            ylab = ylab,
           cex.lab = 1.5,
           cex.axis = 1.2,
           cex.names = 1.2)
         # cex.main = cex)
           # cex.axis = cex)
    # vioplot(df_results, 
    #         col = MODEL_COLORS,  names = MODEL_NAMES,
    #         main = main_title,
    #         xlab = 'Model', 
    #         ylab = 'Posterior Model Probability',
    #         cex.lab = cex,
    #         cex.main = cex, #cex)
    #         cex.axis = cex, 
    #         cex.names = cex, 
    #          cex.sub = cex)
           
    #par(old_par) 
    
  } else {
    par(mar = c(5.5, 5, 4, 4.7)) #rep(4.5, 4)) #c(1.5, 5, 4, 1.5)) #bottom, left, top, right
    par(oma = c(1, 1, 6, 1))
    boxplot(df_results, col = MODEL_COLORS, names = MODEL_NAMES,
            main = main_title, 
            xlab = 'Model', ylab = 'Posterior Model Probability',
            cex.lab = cex, cex.axis = cex, cex.main = cex, 
            cex.names = cex, cex.sub = cex)
  }
  
  # Adjust axis labels size
  #axis(1, cex.axis = cex)
  #axis(2, cex.axis = cex)
  
  if(PDF){
    dev.off()
  }

  #legend('topright', legend = colnames(df_results), fill = box_colors)
  
}

#********************************
#* POSTERIOR PROBABILITY RESULTS
BAR_PLOT_POST_PROBS_PDF <- function(list_vec_results, title = '', RESULTS_FOLDER, data_type, 
                                    PDF = TRUE,
                                    plot_width = 7.0, plot_height = 5.3, cex = 1.0,
                                    MODEL_NAMES = c('Baseline','SSE','SSI', 'SSEB', 'SSIB'),
                                MODEL_COLORS = c('#FFD700', '#6BA6E9',
                                                 '#FF8000', '#6AA84F', '#DC143C')) {
  
  #PLOT MODEL PROBABILITIES
  #PDF
  if(PDF){
    plot_folder = paste0(RESULTS_FOLDER, 'plots/')
    create_folder(plot_folder)
    time_stamp = GET_CURRENT_TIME_STAMP()
    pdf_file = paste0('post_prob_', data_type,  '_mcmc_results_', '_', time_stamp, '.pdf')  
    pdf(paste0(plot_folder, pdf_file), width = plot_width, height = plot_height) 
    
  }
  
  #SETUP
  par(mar = c(4.5,5,4,4))
  par(oma = c(1, 1, 5, 1)) 
  
  #SETUP
  title = bquote(paste('Posterior Model Probabilities - ',.(title)))
  #title = bquote(paste(.(title),' Posterior Model Probabilities. '))
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  
  barplot(t(df_results), beside = TRUE, col = MODEL_COLORS, ylim = c(0, 1),
          names.arg = MODEL_NAMES, 
          main = title,
          xlab = 'Model', ylab = 'Posterior Probability',
          cex.lab=cex+0.25, cex.axis=cex, 
          cex.main= cex + 0.2, 
          cex.names = cex + 0.2, cex.sub=cex)
  #legend('topright', legend = colnames(df_results), fill = MODEL_COLORS) #, cex = 1.1)
  
  dev.off()
}


BAR_PLOT_POST_PROBS <- function(list_vec_results, title = '', cex = 1.4,
                                MODEL_NAMES = c('Base','SSE','SSI', 'SSEB', 'SSIB'),
                                MODEL_COLORS = c('#FFD700', '#6BA6E9',
                                                  '#FF8000', '#6AA84F', '#DC143C')) {
  
  #PLOT MODEL PROBABILITIES
  
  #SETUP
  title = bquote(paste(.(title),' Posterior Model Probabilities '))
  df_results <- as.data.frame(do.call(cbind, list_vec_results))
  
  barplot(t(df_results), beside = TRUE, 
          col = MODEL_COLORS, ylim = c(0, 1),
          names.arg = MODEL_NAMES, main = title,
          xlab = 'Model', ylab = 'Posterior Probability',
          cex.lab=cex+0.2, cex.axis=cex+0.2, cex.main=cex+0.3, cex.sub=cex+0.2,
          cex.names = cex + 0.2)

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


#*******************88
#FUNCTION
BOX_PLOT_MC_ERROR <- function(df_mce_column, model_num, ylim){
  
  #MODELS
  num_runs = length(df_mce_column)
  MODEL_NAMES = c('Baseline', 'SSE', 'SSI', 'SSEB', 'SSIB')
  MODEL_COLORS <- c('#FFD700', '#6BA6E9', '#FF8000', '#6AA84F', 'red')
  
  boxplot(df_mce_column, col = MODEL_COLORS[model_num], #names = model_names,
          main = paste0(MODEL_NAMES[model_num], ' Model - Monte Carlo Error. N = ', num_runs ,
                        ', sd = ', round(sd(df_mce_column), 2)),
          xlab =  paste0(MODEL_NAMES[model_num], ' Model'), ylab = 'log Model Evidence',
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, 
          cex.names = 1.3, cex.sub = 1.3,
          ylim = ylim) 
  
}
