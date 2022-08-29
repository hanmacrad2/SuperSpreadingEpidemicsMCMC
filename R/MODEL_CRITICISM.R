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
#' 
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

#'Calculate a set of summary statistics of the data
#' Returns the summary statistics of the data in a dataframe or list format
#' 
#' @param data Data for which the summary statistics are calculated
#' @param FLAG_MAKE_DF Flag to indicate whether to return a dataframe or a list of summary statistics 
#' @return summary statistics of the data in a dataframe or list format 
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' summary_stats_results = GET_SUMMARY_STATS(data, FLAG_MAKE_DF = TRUE)
#' 
GET_SUMMARY_STATS <- function(data, FLAG_MAKE_DF){

  #Initialise
  mid_point = length(data)/2; end = length(data)
  
  if (FLAG_MAKE_DF){
    
    #Data-frame of summary statistics (20)
    summary_stats_results = data.frame( 
      
      sum_infects = sum(data),
      sum_1st_half  = sum(data[1:mid_point]),
      sum_2nd_half =  sum(data[mid_point+1:end]),
      
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
      
    )
    
  } else {
    
    #List if Dataframe already created
    summary_stats_results = list(sum(data), sum(data[1:mid_point]), sum(data[mid_point +1:end]),
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

####################
RUN_MCMC_REPS <- function(epidemic_data, root_folder,
                          model_type = list(FLAG_BASE = TRUE, FLAG_SSE = FALSE, FLAG_SSI = FALSE),
                          modelling_specs = list(n_reps = 500, n_mcmc = 500000)){
  
  'Run mcmc for n_reps iterations and save' 
  
  for(rep in 1:modelling_specs$n_reps) {
    
    #Folder 
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(base_folder, '/rep_', rep)
    ifelse(!dir.exists(file.path(folder_rep)), dir.create(file.path(folder_rep), recursive = TRUE), FALSE)
    
    #MCMC 
    if (model_type$FLAG_BASE){
      
      mcmc_output = 
    } else if (model_type$FLAG_SSE) {
      
      mcmc_output = SSE_MCMC_ADAPTIVE(epidemic_data)
      
    } else if (model_type$FLAG_SSI) {
      
      #Data (Fix to SSI?)
      mcmc_output = SSI_MCMC_ADAPTIVE(epidemic_data)
    }
    
    mcmc_params = mcmc_sse(sim_data, n, sigma, thinning_factor, folder_rep, rep, burn_in)
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_params, file = paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
  }
  
  
}

run_mcmc_reps_ss <- function(n, n_reps, model_params, sigma, flag_dt, base_folder, burn_in){
  
  'Run mcmc for n_reps iterations and save'
  
  #Get model parameters
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  
  #Data_type
  flag1 = flag_dt[1]; flag2 = flag_dt[2]; flag3 = flag_dt[3] 
  start_rep = 1
  cat('r0 = ', r0, '\n');
  
  # if (flag1){
  #   start_rep = 462
  # } else {
  #   start_rep = 1
  # }
  
  #Repeat for n reps
  for(rep in start_rep:n_reps) {
    
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(base_folder, '/rep_', rep)
    ifelse(!dir.exists(file.path(folder_rep)), dir.create(file.path(folder_rep), recursive = TRUE), FALSE)
    
    #MCMC 
    mcmc_params = mcmc_sse(sim_data, n, sigma, thinning_factor, folder_rep, rep, burn_in)
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_params, file = paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
  }
}

##############################
#GET_SUM_STATS_TOTAL <- function()
  
get_sum_stats_total <- function(base_folder_current, thinning_factor, n_reps, n_mcmc){
  
  'Get summary stats and p vals for all mcmc reps'
  
  for(rep in 1:n_reps) {
    
    #Get results
    folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
    cat('rep = ', rep)
    true_rep_sim = readRDS(paste0(folder_rep, '/sim_data.rds'))
    mcmc_params <- readRDS(paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
    #Get true summary statistics 
    df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    #Save 
    saveRDS(df_true_ss, file = paste0(folder_rep, 'df_true_sum_stats_rep_', rep, '.rds' ))
    
    #Get parameters
    alpha_mcmc = mcmc_params[1]; alpha_mcmc = unlist(alpha_mcmc)
    beta_mcmc = mcmc_params[2]; beta_mcmc = unlist(beta_mcmc)
    gamma_mcmc = mcmc_params[3]; gamma_mcmc = unlist(gamma_mcmc)
    r0_mcmc = mcmc_params[4]; r0_mcmc = unlist(r0_mcmc)
    
    #Simulate data using thinned params
    for(i in seq(burn_in, n_mcmc, by = thinning_factor)){
      #print(paste0("mcmc summary stat rep ", i))
      #Simulate data
      sim_data_model_crit = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alpha_mcmc[i], beta_mcmc[i], gamma_mcmc[i])
      #Save data
      saveRDS(sim_data_model_crit, file = paste0(folder_rep, 'mcmc/sim_data_iter_', i, '.rds' ))
      
      #Get summary stats. 
      if (i == burn_in) { #first rep
        #cat('CREATE  df_summary_stats')
        flag_create = TRUE
        df_summary_stats = get_summary_stats(sim_data_model_crit, flag_create)
        flag_create = FALSE 
        #Get indices of iterations
        list_ss_iters = c(i)
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats(sim_data_model_crit, flag_create)
        list_ss_iters = c(list_ss_iters, i)
      }
      
    }
    
    #Save summary stats
    saveRDS(df_summary_stats, file = paste0(folder_rep, '/df_summary_stats_', rep, ".rds"))
    #print(paste0('df_summary_stats', df_summary_stats))
    #Save ss iterations
    saveRDS(list_ss_iters, file = paste0(folder_rep, '/list_ss_iters_i', rep, '.rds'))  
    
  }
}


#2B. TOTAL SUMMARY STATS 
get_sum_stats_total <- function(base_folder_current, thinning_factor, n_reps, n_mcmc){
  
##############################
#3. GET P VALUES FOR ALL  REPS
get_p_values_total <- function(base_folder_current, n_reps){ }