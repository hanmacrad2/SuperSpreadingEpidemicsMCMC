#RUN MCMC CREDIBLE INTERVALS
source("R/UTIL_FUNCTIONS.R")
source("R/SSID_MCMC_ADAPTIVE.R")

#OUTPUT
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/ssid_model/v1"
  
#PARAMS
num_days = 110
alpha_vec = seq(from = 0.9, to = 2, length = 8)

#MCMC FOR EACH VALUE
par(mfrow = c(4,5))
for(i in seq_along(alpha_vec)){
  print(paste0('i: ', i))
  
  #1. SIMULATE DATA
  print(alpha_vec[i]);
  r0 = alpha_vec[i] + 0.5
  data_sse = SIMULATION_SSE(alpha_vec[i])
  plot.ts(data_sse, main = print(paste0(round(r0, 2))))
}

#MCMC SSE POISSON COMPOUND
seedX = 4
set.seed(4)
data = SIMULATE_SSID(); dataI = data$epidemic_data
plot.ts(dataI)

#START MCMC
Rprof(tmp <- tempfile())
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssid = MCMC_ADAPTIVE_SSID(dataI, OUTER_FOLDER, seedX)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssid$time_elap = time_elap
Rprof()
summaryRprof(tmp)

#PLOT
df = PLOT_SSID_MCMC_GRID(dataI, mcmc_ssid, seedX, data$eta_vec, -500)

#*********************************************
#SETUP!
alpha_vec = seq(from = 0.9, to = 2, length = 8); 
num_iters = length(alpha_vec)
vec_means = vector("numeric", length = num_iters)
vec_lower = vector("numeric", length = num_iters); vec_upper = vector("numeric", length = num_iters)

vec_meansb = vector("numeric", length = num_iters)
vec_lowerb = vector("numeric", length = num_iters); vec_upperb = vector("numeric", length = num_iters)

#RUN!
for(i in 1:num_iters){
  print(paste0('i: ', i))
  
  #1. SIMULATE DATA
  set.seed(i); print(alpha_vec[i]);
  data_ssidX = SIMULATE_SSID(alphaX = alpha_vec[i])
  data_ssid = data_ssidX$epidemic_data 
  plot.ts(data_ssid)
  saveRDS(data_ssid, file = paste0(OUTER_FOLDER, '/data_ssid_', i, '.rds' ))
  #SAVE DATA
  print('1')
  #2. MCMC 
  start_time = Sys.time()
  print(paste0('start_time:', start_time))
  mcmc_ssid_output = MCMC_ADAPTIVE_SSID(data_ssid, OUTER_FOLDER, i)
  end_time = Sys.time()
  time_elap = get_time(start_time, end_time)
  mcmc_ssid_output$time_elap = time_elap
  saveRDS(mcmc_ssid_output, file = paste0(OUTER_FOLDER, '/mcmc_ssid_', i, '.rds' ))
  print('2')
  #MEANS
  vec_means[seedX] = mean(mcmc_ssid_output$ssid_params_matrix[, 1])
  #vec_lower[seedX] = get_lower_ci(mcmc_ssid_output$ssid_params_matrix[, 1]) # get_lower_ci(mcmc_ssid_output$ssid_params_matrix[, 1]) 
  #vec_upper[seedX] = get_upper_ci(mcmc_ssid_output$ssid_params_matrix[, 1]) 
  
  #vec_meansb[seedX] = mean(mcmc_ssid_output$ssid_params_matrix[, 2])
  #vec_lowerb[seedX] = get_lower_ci(mcmc_ssid_output$ssid_params_matrix[, 2]) 
  vec_upperb[seedX] = get_upper_ci(mcmc_ssid_output$ssid_params_matrix[, 2]) 
  
  #PRINT
  print('alpha')
  print(vec_means); print(vec_lower); print(vec_upper)
  print('k')
  print(vec_meansb); print(vec_lowerb); print(vec_upperb)
}
