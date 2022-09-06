#DISPLAY MODEL CRITICISM RESULTS

#***************
#PLOT P VALUES
#' Plot histograms of the ppp-values for each summary statistic 
#'
#' Plot histograms of the ppp-values (posterior predictive p-values) for each of the summary statistics collected in the model criticism procedure
#'
#' @param df_p_values dataframe of the ppp-values
#' @param model_type Type of epidemic model fitted to the data
#' 
#' @return Plots of histograms of the ppp-values  
#' @export PLOT_P_VALUES
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples PLOT_P_VALUES(df_p_values, 'SSE')

PLOT_P_VALUES <- function(FOLDER_REP, model_type){ #*OR ROOT FOLDER
  
  'Plot histograms of the ppp-values'
  df_p_vals =  readRDS(paste0(FOLDER_REP, '/df_total_p_values.rds'))
  
  par(mfrow=c(4,5)) #c(3,4)
  num_iters = length(df_p_values[,1])
  val_05 = 0.05
  
  for (i in c(1:20)){ 
    
    #ppp-value <0.05
    percent_lt_05 = (length(which(df_p_values[,i] < val_05))/num_iters)*100
    
    #Histogram
    if (i == 1){
      print(paste0('i = ', i))
      hist(df_p_vals[,i], breaks = 100, #freq = FALSE, 
           xlab = paste0('p value, < 0.05: ', percent_lt_05, '%'),
           ylab = 'Num Samples', col = 'green',
           main = paste('', toupper(colnames(df_p_vals)[i]),', Model:', model_type),
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(v = 0.05, col = 'red', lwd = 2)
      
    } else {

    hist(df_p_values[,i], breaks = 100, #freq = FALSE, 
         xlab = paste0('p value, < 0.05: ', percent_lt_05, '%'),
         ylab = 'N Samples', col = 'green',
         main = paste('', toupper(colnames(df_p_values)[i]), ', n reps:', num_iters),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = 0.05, col = 'red', lwd = 2)
    
    }
  }
}

#***************
#PLOT SUMMARY STATS
#' Plot histograms of the ppp-values for each summary statistic 
#'
#' Plot histograms of the ppp-values (posterior predictive p-values) for each of the summary statistics collected in the model criticism procedure
#'
#' @param df_p_values dataframe of the ppp-values
#' @param model_type Type of epidemic model fitted to the data
#' 
#' @return Plots of histograms of the ppp-values  
#' @export PLOT_P_VALUES
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples PLOT_SUMMARY_STATS(df_p_values, 'SSE')
# 
# PLOT_SUMMARY_STATS <- function(true_r0, model_type, sim_data_rep, df_summary_stats, df_true_sum_stats,
#                                list_p_vals, upper_quant, trim_flag){
  
PLOT_SUMMARY_STATS <- function(FOLDER_REP, epidemic_data, data_type, rep) {
'Plot sim data, summary stats and true summary stat for a given mcmc rep' 
  
  #DATA (Match to TERMINAL!) 
  list_p_vals <- readRDS(paste0(FOLDER_REP, 'list_p_vals_rep', rep, '.rds'))
  df_summary_stats <- readRDS(paste0(FOLDER_REP, 'df_replicated_summary_stats_', rep, '.rds')) 
  df_true_sum_stats <- readRDS(paste0(FOLDER_REP, 'df_ss_true_rep_', rep, '.rds')) 
  len_data = length(list_p_vals)
  
  #PLOT SETUPT
  par(mfrow = c(5, 5)) 
  colorsX <- rainbow(len_data+1); colors_line <- rainbow(c(15:15+len_data+1))
  
  #PLOTS
  #i. EPIDEMIC DATA
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count', 
          main = paste0(data_type, "epidemic data"), 
          col = colorsX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. SUMMARY STATISTICS
  for (i in 1:len_data){
    
    X = df_summary_stats[,i] 
    
    if (!is.na(X)) {
      
      #HISTOGRAM
      hist(df_summary_stats[,i], breaks = 100, #freq = FALSE, 
           main = paste('', toupper(colnames(df_summary_stats)[i]),', p value:', round(list_p_vals[i],3)),
           xlab = paste0('', toupper(colnames(df_summary_stats)[i]), ', T:',
                         round(df_true_sum_stats[nrow(df_true_sum_stats), i], 2)),
           ylab = 'Num Samples',
           col = colorsX[i+1],
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(v = df_true_sum_stats[1, i], col = colors_line[i], lwd = 2.5) #This should be true value, not p value
    }

    
    # if (trim_flag){
    #   print('trimmed')
    #   X = upper_quantile(X, upper_quant)
    # }
    
    # print(paste0('Col i = ', colnames(df_summary_stats)[i]))
    # print(paste0('True val i = ', round(df_true_sum_stats[nrow(df_true_sum_stats), i], 2)))
    # print(paste0('X = ', X))
    
  }
  
}

#'
#'
#'
PLOT_REPLICATED_DATA <- function(FOLDER_REP, epidemic_data, rep, data_type){
  
  #DATA REPLICATED
  FOLDER_REP_DATA = paste0(FOLDER_REP, 'replicated_data/')
  print(FOLDER_REP_DATA)
  mcmc_output <- readRDS(paste0(FOLDER_REP, 'mcmc_output_rep_', rep, '.rds'))
  list_rep_idx <- readRDS(paste0(FOLDER_REP, 'list_ss_iters_i', rep, '.rds'))
  sample_reps =  sort(sample(list_rep_idx, size = 11)); coloursX <- rainbow(length(sample_reps)+1)
 
  #PLOTS
  par(mfrow=c(3,4))
  #i. EPIDEMIC DATA
  plot.ts(epidemic_data, xlab = 'Time', ylab = 'Daily Infections count', 
          main = paste0(data_type, " Epidemic data"), 
          col = coloursX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. REPLICATED DATA
  for (i in 1:length(sample_reps)){
    
    print(paste0(FOLDER_REP_DATA, 'sim_data_iter_', sample_reps[i], '.rds'))
    #DATA
    replicated_data = readRDS(paste0(FOLDER_REP_DATA, 'sim_data_iter_', sample_reps[i], '.rds'))
    #PLOT
    plot.ts(replicated_data, xlab = 'Time', ylab = 'Daily Infections count',
            main = paste0("Iter_", i, 'R0:', round(mcmc_output$r0_vec[i], 3)),
            col = coloursX[i], lwd = 2,
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  
}