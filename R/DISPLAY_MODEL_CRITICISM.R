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

PLOT_P_VALUES <- function(df_p_values, model_type){
  
  'Plot histograms of the ppp-values'
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
#' @examples PLOT_P_VALUES(df_p_values, 'SSE')

plot_summary_stats <- function(true_r0, model_type, sim_data_rep, df_sum_stats, df_true_sum_stats,
                               list_p_vals, upper_quant, trim_flag){
  
  'Plot sim data, summary stats and true summary stat for a given mcmc rep' 
  #Setup
  par(mfrow = c(5, 5)) 
  len_data = length(list_p_vals)
  print(paste0('length of data:', len_data))
  colorsX <- rainbow(len_data+1)
  colors_line <- rainbow(c(15:15+len_data+1))
  
  #i. Sim_data - Infections
  plot.ts(sim_data_rep, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(rep, "Infts,", model_type, "R0 = ", true_r0), #model_type
          col = colorsX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Columns
  for (i in c(1:len_data)){
    print(paste0('i = ', i))
    X = df_sum_stats[,i] # df_sum_stats[1:nrow(df_sum_stats),i]
    
    if (trim_flag){
      print('trimmed')
      X = upper_quantile(X, upper_quant)
    }
    
    # print(paste0('Col i = ', colnames(df_sum_stats)[i]))
    # print(paste0('True val i = ', round(df_true_sum_stats[nrow(df_true_sum_stats), i], 2)))
    # print(paste0('X = ', X))
    
    #Histogram
    hist(X, breaks = 100, #freq = FALSE, 
         xlab = paste0('', toupper(colnames(df_sum_stats)[i]), ', T:', round(df_true_sum_stats[nrow(df_true_sum_stats), i], 2)),
         ylab = 'Num Samples',
         col = colorsX[i+1],
         main = paste('', toupper(colnames(df_sum_stats)[i]),', p value:', round(list_p_vals[i],3)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = df_true_sum_stats[nrow(df_true_sum_stats), i], col = colors_line[i], lwd = 2.5) #This should be true value, not p value
    
  }
  
}