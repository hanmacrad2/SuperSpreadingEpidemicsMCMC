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
#' @param FLAG_DF_CREATE Flag to indicate whether to return a dataframe or a list of summary statistics 
#' @return summary statistics of the data in a dataframe or list format 
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#' summary_stats_results = GET_SUMMARY_STATS(data, FLAG_DF_CREATE = TRUE)
#' 
GET_SUMMARY_STATS <- function(data, FLAG_DF_CREATE = TRUE){

  #Initialise
  mid_point = length(data)/2; end = length(data)
  
  if (FLAG_DF_CREATE){
    
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

#' Run MCMC sampler for a number of reps to facilitate Model Criticism 
#'
#' MCMC algorithm with Adaptation for obtaining samples from the parameters of a
#' super-spreading events (SSE) epidemic model
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param root_folder root folder location in which to store the mcmc results for \code{"n_reps"}
#' @param model_type A list of indicators of type boolean to indicate the type of epidemic model to fit and corresponding MCMC sampler to run;
#' \itemize{
#'   \item \code{"FLAG_BASE"} - if TRUE the baseline model MCMC sampler is implemented 
#'   \item \code{"FLAG_SSE"} - if TRUE the super-spreading events (SSE) model MCMC sampler is implemented 
#'   \item \code{"FLAG_SSI"} - if TRUE the super-spreading individuals (SSI) model MCMC sampler is implemented 
#' @param modelling_specs A list of specifications;
#' \itemize{
#'   \item \code{"n_reps"} - Number of repetitions of the mcmc sampler. Multiple runs required for model criticism
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler }
#' @return saves the \code{"mcmc_output"} int the \code{"root_foler"}/\code{"rep"} location
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#'RUN_MCMC_REPS(epidemic_data, results_folder)
#' 
RUN_MCMC_REPS <- function(epidemic_data, root_folder,
                          model_type = list(FLAG_BASE = TRUE, FLAG_SSE = FALSE, FLAG_SSI = FALSE),
                          modelling_specs = list(n_reps = 500, n_mcmc = 500000, burn_in_factor = )){
  
  'Run mcmc for n_reps iterations and save' 
  
  for(rep in 1:modelling_specs$n_reps) {
    
    #Folder 
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(base_folder, '/rep_', rep)
    ifelse(!dir.exists(file.path(folder_rep)), dir.create(file.path(folder_rep), recursive = TRUE), FALSE)
    
    #MCMC 
    if (model_type$FLAG_BASE) {
      mcmc_output = BASELINE_MCMC_ADAPTIVE(epidemic_data)
      
    } else if (model_type$FLAG_SSE) {
      mcmc_output = SSE_MCMC_ADAPTIVE(epidemic_data)
      
    } else if (model_type$FLAG_SSI) {
      #Data (Fix to SSI?)
      mcmc_output = SSI_MCMC_ADAPTIVE(epidemic_data)
    }
    
    mcmc_output = mcmc_sse(sim_data, n, sigma, thinning_factor, folder_rep, rep, burn_in)
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_output, file = paste0(folder_rep, '/mcmc_output_rep_', rep, '.rds' ))
  }
}

#GET_SUM_STATS_TOTAL
GET_SUM_STATS_TOTAL <- function(epidemic_data, root_folder,
                                 model_type = list(FLAG_BASE = TRUE, FLAG_SSE = FALSE, FLAG_SSI = FALSE),
                                 modelling_specs = list(n_reps = 500, n_mcmc = 500000, thinning_factor = 10, 
                                                        burn_in_size = 0.01)) { #thinning factor?
  
  'Get summary stats and p valuess for all mcmc reps'

  #MODELLING SPECS
  burn_in = modelling_specs$burn_in_size*modelling_specs$n_mcmc
  num_days = length(epidemic_data)
  
  for(rep in 1:modelling_specs$n_reps) {
    
    #GET RESULTS
    print(paste0('rep = ', rep))
    folder_rep = paste0(root_folder, "/rep_", rep, '/')
    
    #MCMC
    mcmc_output <- readRDS(paste0(folder_rep, '/mcmc_output_rep_', rep, '.rds' ))
    
    #GET SUMMARY STATS (TRUE)
    df_ss_true = GET_SUMMARY_STATS(epidemic_data) 
    saveRDS(df_ss_true, file = paste0(folder_rep, 'df_ss_true_rep_', rep, '.rds' )) #save
    
    #GET MCMC PARAMETERS
    # if (model_type$FLAG_BASE) {
    #   mcmc_output = BASELINE_MCMC_ADAPTIVE(epidemic_data)
    # 
    # } else if (model_type$FLAG_SSE) {
    #   mcmc_output = SSE_MCMC_ADAPTIVE(epidemic_data)
    # 
    # } else if (model_type$FLAG_SSI) {
    #   #Data (Fix to SSI?)
    #   mcmc_output = SSI_MCMC_ADAPTIVE(epidemic_data)
    # }
    # 
    # #Get parameters -> in function
    # alpha_mcmc = mcmc_output[1]; alpha_mcmc = unlist(alpha_mcmc)
    # beta_mcmc = mcmc_output[2]; beta_mcmc = unlist(beta_mcmc)
    # gamma_mcmc = mcmc_output[3]; gamma_mcmc = unlist(gamma_mcmc)
    # r0_mcmc = mcmc_output[4]; r0_mcmc = unlist(r0_mcmc)
    
    #SIMULATE DATA (THINNED)
    for(i in seq(burn_in, modelling_specs$n_mcmc, by = modelling_specs$thinning_factor)){
      
      if (model_type$FLAG_BASE) {
        R0 = mcmc_output$r0_vec[i]
        sim_data = SIMULATE_BASELINE_EPIDEMIC(R0, num_days = num_days) 
      }  else if (model_type$FLAG_SSE) {
        
        sim_data = SIMULATE_SS_EVENTS(mcmc_output$alpha_vec[i], mcmc_output$beta_vec[i], mcmc_output$gamma_vec[i])

      } else if (model_type$FLAG_SSI) {

        sim_data = SSI_MCMC_ADAPTIVE(mcmc_output$a_vec[i], mcmc_output$b_vec[i], mcmc_output$c_vec[i])
      }

      #Save data
      saveRDS(sim_data, file = paste0(folder_rep, 'mcmc/sim_data_iter_', i, '.rds' ))
      
      #GET SUMMARY STAT RESULTS
      if (i == burn_in) { #first rep
        df_summary_stats = GET_SUMMARY_STATS(sim_data)
        list_ss_iters = c(i)
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = GET_SUMMARY_STATS(sim_data, FLAG_DF_CREATE = FALSE)
        list_ss_iters = c(list_ss_iters, i)
      }
      
    }
    
    #SAVE SUMMARY STATISTICS + ITERATIONS
    saveRDS(df_summary_stats, file = paste0(folder_rep, '/df_summary_stats_', rep, ".rds"))
    saveRDS(list_ss_iters, file = paste0(folder_rep, '/list_ss_iters_i', rep, '.rds'))  
  }
}

#
#3. GET P VALUES FOR ALL  REPS
get_p_values_total <- function(root_folder, n_reps){ }

#####
get_p_values_total <- function(base_folder_current, n_reps){
  
  for(rep in 1:n_reps) {
    cat('rep = ', rep)
    #Get results
    folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
    cat('folder_rep', folder_rep)
    true_rep_sim = readRDS(paste0(folder_rep, '/sim_data.rds'))
    #Get true summary statistics 
    df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    
    #Data
    df_summary_stats_rep <- readRDS(paste0(folder_rep, '/df_summary_stats_', rep, '.rds' ))
    
    #Get p values
    list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_p_vals, file = paste0(folder_rep, '/list_p_vals_rep', rep, ".rds"))
    
    list_all_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values_list(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_all_p_vals, file = paste0(folder_rep, '/list_all_p_vals_rep_', rep, ".rds"))
    
    #Save all 
    if (!exists("df_p_values")) {
      df_p_values = data.frame(sum_infects = list_p_vals[1],
                               sum_1st_half = list_p_vals[2],
                               sum_2nd_half = list_p_vals[3],
                               median_infect = list_p_vals[4],
                               max_infect = list_p_vals[5],
                               sd_infect = list_p_vals[6],
                               infect_q_75 = list_p_vals[7],
                               infect_q_87_5 = list_p_vals[8],
                               max_dif = list_p_vals[9],
                               med_dif  = list_p_vals[10],
                               max_dif2nd =  list_p_vals[11],
                               med_dif2nd =  list_p_vals[12],
                               max_dif_normI =  list_p_vals[13],
                               max_dif_normII =  list_p_vals[14],
                               max_dif_normIII = list_p_vals[15],
                               max_dif2nd_I =  list_p_vals[16],
                               max_dif_2ndII =  list_p_vals[17],
                               med_dif_normI =  list_p_vals[18],
                               med_dif_normII =  list_p_vals[19],
                               med_dif_normIII =  list_p_vals[20]
                               
      )
      print(paste0('df_p_values', df_p_values))
      
    } else {
      df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
    }
    
  }
  
  #Ensure its a df
  df_p_values = as.data.frame(df_p_values)
  
  #SaveRDS
  saveRDS(df_p_values, file = paste0(base_folder_current, '/total_p_values_iter_', iter, '.rds' ))
  
  #Return p values
  df_p_values
  
}