#SSEB MODEL + REAL DATA

DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'
filename = 'df_nz_south_20.rds'
df_epi_data = readRDS(paste0(DATA_FOLDER, filename))

epidemic_data = df_epi_data$cases
plot.ts(epidemic_data)

MCMC_INFER_SSEB(epidemic_data, 10000)


#CHECK
num_days = length(epidemic_data)

t = 2
shape_gamma = 6; scale_gamma = 1
r0 = 2; alpha = 0.3; beta = 8.6
gamma = r0*(1 - alpha)/beta
prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
  pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)

for (t in 2:num_days) {
  
    xt = epidemic_data[t]
    #PROBABILITY_XT(epidemic_data[t], lambda_vec[t], gamma, beta, alpha, r0)
  
    max_et = GET_MAX_SS_EVENTS(xt) #ADDED 12/01/24
    prob_xt = 0
    vec_prob_xt = vector('numeric', length = max_et + 1)
    lambda_t = sum(epidemic_data[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    
    gamma_lambda_t = gamma*lambda_t
    alpha_r0_lambda_t = alpha*r0*lambda_t
    
    et_values <- 0:max_et
    vec_prob_et = dpois(et_values, gamma_lambda_t) #, log = TRUE) 
  
}


for (t in 2:num_days) {
  
  #print(t)
  #print(epidemic_data[t])
  lambda_t = sum(epidemic_data[1:(t-1)]*rev(prob_infect[1:(t-1)]))
  prob_xt = PROBABILITY_XT_NON_LOG(epidemic_data[t], lambda_t, gamma, beta, alpha, r0)
  print(round(prob_xt, 4))
  
}

PROBABILITY_XT_NON_LOG <- function(xt, lambda_t, gamma, beta, alpha, r0, max_et = 100) {
  
  'Compound Poisson Prob for SSEB model'
  
  #MAX EVENTS (ADD)
  max_et = GET_MAX_SS_EVENTS(xt) #ADDED 12/01/24
  prob_xt = 0
  vec_prob_xt = vector('numeric', length = max_et + 1)
  
  gamma_lambda_t = gamma*lambda_t
  alpha_r0_lambda_t = alpha*r0*lambda_t
  
  et_values <- 0:max_et
  vec_prob_xt = dpois(et_values, gamma_lambda_t)* dpois(xt, beta*et_values + alpha_r0_lambda_t)
  
  #prob_xt = LOG_SUM_EXP(vec_prob_xt)
  
  return(prob_xt)
}

#data_type = 'nz_south_wedding'
file_data = 'epi_data_nz.rds'
epidemic_data = readRDS(file_data)

#data_type = 'nz_south_03_20'
#file_data = 'epi_data_nz_south_03_20.rds'