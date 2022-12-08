#RUN SSEB MCMC 
#FOLDERS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R'); source('R/EPI_FUNCTIONS.R')
source('R/SSEB_MCMC.R')

#OUTPUT FOLDER
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/sseb"
seedX = 10

#CREATE OUTPUT FOLDER
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', seedX)
ifelse(!dir.exists(file.path(CURRENT_OUTPUT_FOLDER)),
       dir.create(file.path(CURRENT_OUTPUT_FOLDER), recursive = TRUE), FALSE)

#PARAMETERS 
n_mcmc = 50000 
alphaX = 0.8; num_days = 50 #110 #Question: 50 days ok? Run 110 now. 1000 runs takes x mins

#1.SIMULATE DATA
set.seed(seedX)
epi_data_sseb = SIMULATE_EPI_SSEB(num_days)
plot.ts(epi_data_sseb)
#SAVE DATA
saveRDS(epi_data_sseb, file = paste0(CURRENT_OUTPUT_FOLDER, '/epi_data_sseb_', seedX, '.rds' ))
epi_data_sseb = sim10
#2.RUN MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sseb_output10 = MCMC_INFER_SSEB(epi_data_sseb, n_mcmc)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sseb_output10$time_elap = time_elap

#SAVE MCMC RESULTS
saveRDS(mcmc_sseb_output10, file = paste0(CURRENT_OUTPUT_FOLDER, '/mcmc_sseb_', seedX, '.rds' ))

#SEED
seedX = seedX + 1

########################
#CREDIBLE INTERVALS
alpha_vec 
alpha_vec = c(0.7, 0.8, 0.9, 1.0, 1.1, 0.7, 0.8, 0.9, 1.0, 1.1)
beta_vec = c(.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1)
r0_vec = alpha_vec + 10.*beta_vec

i = 1

i = i + 1
i
#Simulate data
set.seed(i)
sim10 =  SIMULATE_EPI_SSEB(num_days, alphaX = alpha_vec[i], betaX = beta_vec[i])
plot.ts(sim10)

#Data
list_epi_data = list(sim1 = sim1, sim2 = sim2, sim3 = sim3,
                     sim4 = sim4, sim5 = sim5,
                     sim6 = sim6, sim7 = sim7, sim8 = sim8,
                     sim9 = sim9, sim10 = sim10)

#PLOT
par(mfrow = c(3,4))
for(i in seq_along(alpha_vec)){
  print(paste0('i: ', i))
  plot.ts(list_epi_data[i][1],
          main = paste0('R0: ', r0_vec[i], 'alpha: ', alpha_vec[i], 'beta: ', beta_vec[i]))
}


epi_data_sseb = SIMULATE_EPI_SSEB(num_days, alphaX = alpha_vec[i], alphaX = alpha_vec[i],)

for(i in seq_along(alpha_vec)){
  print(paste0('i: ', i))
  
  #1. SIMULATE DATA
  set.seed(i); print(alpha_vec[i]);
  data_sse = SIMULATION_SSE(alpha_vec[i])
  #SAVE DATA
  saveRDS(data_sse)
  
  #2. MCMC
  start_time = Sys.time()
  print(paste0('start_time:', start_time))
  mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(data_sse)
  end_time = Sys.time()
  time_elap = get_time(start_time, end_time)
  mcmc_sse_output$time_elap = time_elap
  saveRDS(mcmc_sse_output)
  
  #MEANS
  vec_means[seedX] = mean(mcmc111$nu_params_matrix[, 1])
  vec_lower[seedX] = get_lower_ci(mcmc111) 
  vec_upper[seedX] = get_upper_ci(mcmc111) 
  
  #PRINT
  print(vec_means); print(vec_lower); print(vec_upper)
}