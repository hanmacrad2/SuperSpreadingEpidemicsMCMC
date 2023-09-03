#RUN MCMC CREDIBLE INTERVALS
setwd('~/GitHub/SuperSpreadingEpidemicsMCMC/')
source('R/UTIL_FUNCTIONS.R')
source('R/SSE_MCMC_POISSON_COMPOUND.R')

#OUTPUT
OUTER_FOLDER = ""
  
#PARAMETERS
num_days = 110
alpha_vec = seq(from = 0.7, to = 3, length = 50) #TIME OF 1000
alphaX = 1.4
data_sse = SIMULATION_SSE(alpha_vec[i])


#APPLY MULTIPLE TIME - CREDIBLE INTERVALS 
#MCMC FOR EACH VALUE
par(mfrow = c(4,5))

for(i in seq_along(alpha_vec)){
  print(paste0('i: ', i))
  
  #1. SIMULATE DATA
  print(alpha_vec[i]);
  r0 = alpha_vec[i] + 0.5
  data_sse = SIMULATION_SSE(alpha_vec[i])
  #SAVE DATA
  plot.ts(data_sse, main = print(paste0(round(r0, 2))))
  
}

#SETUP!
num_iters = length(alpha)
vec_means = vector("numeric", length = num_iters)
vec_lower = vector("numeric", length = num_iters); vec_upper = vector("numeric", length = num_iters)


#RUN!
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



#DATA
seedX = 1; set.seed(seedX)
data_sse = SIMULATION_SSE(1.0)
plot.ts(data_sse)

#MCMC SSE POISSON COMPOUND
#START MCMC
Rprof(tmp <- tempfile())
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_sse_output = SSE_POI_MCMC_ADAPTIVE(data_sse)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_sse_output$time_elap = time_elap
Rprof()
summaryRprof(tmp)

#Plot
#2. SSE
df_sse = PLOT_SS_MCMC_GRID(data_sse, mcmc_sse_output,
                           mcmc_specs = list(model_type = 'SSE', n_mcmc = 100, 
                                             mod_start_points = list(m1 = 0.8, m2 = 0.1, m3 = 10),
                                             mod_par_names = c('alpha', 'beta', 'gamma'),
                                             seed_count = 1, burn_in_pc = 0.05, thinning_factor = 10),
                           priors_list = list(a_prior_exp = c(1, 0), b_prior_ga = c(10, 2/100), b_prior_exp = c(0.1,0), #10, 1/100
                                              c_prior_ga = c(10, 1), c_prior_exp = c(0.1,0)),
                           FLAGS_LIST = list(SSI = FALSE, BURN_IN = FALSE, THIN = FALSE,
                                             DATA_AUG = FALSE, ADAPTIVE = TRUE, MULTI_ALG = FALSE,
                                             PRIOR = TRUE, B_PRIOR_GAMMA = FALSE, C_PRIOR_GAMMA = FALSE))
