#SCRIPT

args <- commandArgs(trailingOnly = TRUE)
array_index <- as.integer(args[1])

# Use the array index for your specific task
# For example, you can use it to set a random seed:
set.seed(array_index)

#Loop and repeat
n_mcmc = 30000
num_reps = 10
list_r0 = c(0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,1.8,1.9,2.0)

for (r0 in list_r0){
  
  for (i in 1:num_reps){
    
    data_baseline = SIMULATE_EPI_BASELINE(list_r0[r0])
    filename <- paste("data_baseline_", r0, "_", i, ".rds", sep = "")
    saveRDS(data_baseline, file = filename)
    
    #MCMC
    mcmc_output = MCMC_INFER_BASELINE(data_baseline, n_mcmc = n_mcmc)
    mcmc_file <- paste("mcmc_baseline_", r0, "_", i, ".rds", sep = "")
    saveRDS(data_baseline, file = mcmc_file)
  }
}
