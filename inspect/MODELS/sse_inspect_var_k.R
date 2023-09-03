#SSE MODEL INSPECT VARIANCE VS k

var_neg_bin <- function(k, mean_val = 1.5){
 
  var = mean_val + (mean_val^2)/k
  return(var)
}


#SETUP
plot_neg_bin_var <- function(num_reps = 5000){
  
  list_k = seq(from = 0.05, to = 10, length = num_reps)
  list_var = vector('numeric', length = num_reps)
    
  for (i in 1:length(list_k)){
    
    list_var[i] = var_neg_bin(list_k[i])
    
  }
  
  #Plot
  #ticks <- seq(0.1, 10, by = 0.5)
  #axis(1, at = ticks, labels = ticks)
  
  plot(list_k, list_var, type = 'p',
       main = 'Variance of Neg Bin vs k. Constant Mean = 1.5', xlab = 'k', ylab = 'Variance')

} 

#PLOT
plot_neg_bin_var()


#TEST
# mean_val = 1.5
# k = 0.5
# var = mean_val + (mean_val^2)/k
# var
