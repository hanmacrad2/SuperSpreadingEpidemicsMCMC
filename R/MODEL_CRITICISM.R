#MODEL CRITICISM

#'P value derived by comparing the true summary statistics of the true epidemic data with that of the summary statistics of the simulated data  
#'
#' Returns the p value 
#' 
#' @param column_true_val Column of the true summary stat values of the epidemic data 
#' @param column_summary_stat Column of summary stats of the simulated data
#' @return p-value 
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' p_value = GET_P_VALUE(column_true_val, column_summary_stat)

GET_P_VALUE <- function(column_true_val, column_summary_stat) {
  'Get p values - comparing  summary stat columns to true value'
  
  #Final val
  num_iters = length(column_sum_stat)# - 1
  #P value
  prop_lt = length(which(column_sum_stat < column_true_val))/num_iters + 0.5*(length(which(column_sum_stat == column_true_val)))/num_iters
  prop_gt = length(which(column_sum_stat > column_true_val))/num_iters + 0.5*(length(which(column_sum_stat == column_true_val)))/num_iters
  pvalue = min(prop_lt, prop_gt); pvalue = 2*pvalue
  
  pvalue
  
}

#'P value derived by comparing the true summary statistics of the true epidemic data with that of the summary statistics of the simulated data  
#'
#' Returns the p value 
#' 
#' @param column_true_val Column of the true summary stat values of the epidemic data 
#' @param column_summary_stat Column of summary stats of the simulated data
#' @return p-value 
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' p_value = GET_P_VALUE(column_true_val, column_summary_stat)

GET_SUMMARY_STATS <- function(data, flag_create){
  
  'Calculate summary statistics of the simulated data'
  #Summary stats params
  start_half2 = (length(data)/2)+1
  stop_half2 = length(data)
  
  if (flag_create){
    
    #Df of summary stats (20)
    summary_stats_results = data.frame( 
      
      sum_infects = sum(data),
      sum_1st_half  = sum(data[1:(length(data)/2)]),
      sum_2nd_half =  sum(data[start_half2:stop_half2]),
      
      median_infect = median(data),
      max_infect = max(data),
      sd_infect = sd(data),
      
      infect_q_75 = quantile(data)[4][1][1],
      infect_q_87_5 = quantile(data, probs = seq(0, 1, 0.125))[8][1][1],
      
      #Differences
      max_dif = max((diff(data))), #Change from absolute difference
      med_dif = median(diff(data)),
      max_dif2nd = max(diff(diff(data))),
      med_dif2nd = median(diff(diff(data))),
      
      #Norm Differences
      max_dif_normI = max(diff(data)/(data[1:length(data)-1]+1)),
      max_dif_normII = max(diff(data)/(rollapply(data, 2, mean, by = 1)+1)), #mean of consecutive counts 
      max_dif_normIII = max(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1)), #ratio of consecutive counts
      
      max_dif2nd_I = max(diff(diff(data))/(data[1:length(data)-1]+1)),
      max_dif_2ndII = max(diff(diff(data))/(rollapply(data, 2, mean, by = 1)+1)),
      
      med_dif_normI = median(diff(data)/(data[1:length(data)-1]+1)),
      med_dif_normII = median(diff(data)/(rollapply(data, 2, mean, by = 1) +1)),
      med_dif_normIII = median(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1))
      
      #med_dif2nd_I = median(diff(diff(data))/(data[1:length(data)-1]+1)),
      #med_dif_2ndII = median(diff(diff(data))/(rollapply(data, 2, mean, by = 1)+1))
      
    )
    
  } else {
    #List if df already created
    summary_stats_results = list(sum(data), sum(data[1:(length(data)/2)]), sum(data[start_half2:stop_half2]),
                                 median(data), max(data), sd(data),
                                 quantile(data)[4][1][1], quantile(data, probs = seq(0, 1, 0.125))[8][1][1],
                                 max(diff(data)), median(diff(data)), max(diff(diff(data))), median(diff(diff(data))),
                                 max(diff(data)/(data[1:length(data)-1]+1)), max(diff(data)/(rollapply(data, 2, mean, by = 1)+1)),
                                 max(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1)),
                                 max(diff(diff(data))/(data[1:length(data)-1]+1)),  max(diff(diff(data))/(rollapply(data, 2, mean, by = 1)+1)),
                                 median(diff(data)/(data[1:length(data)-1]+1)), median(diff(data)/(rollapply(data, 2, mean, by = 1)+1)),
                                 median(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1))
                                 
                                 
    )
  }
  
  summary_stats_results
  
}