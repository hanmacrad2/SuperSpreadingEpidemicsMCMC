#MODEL CRITICISM

#' Perform Model Criticism of an epidemic model via posterior predictive checking 
#'
#' Model Criticism of an epidemic model via posterior predictive checking. Calls a series of four functions to execute the process fully. 
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param root_folder root folder location in which to store results
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
#'RUN_MODEL_CRITICISM(epidemic_data, root_folder, model_type = list(FLAG_BASE = FALSE, FLAG_SSE = FALSE, FLAG_SSI = TRUE))
#'
#'
RUN_MODEL_CRITICISM <- function(epidemic_data, root_folder,
                                model_type = list(FLAG_BASE = TRUE, FLAG_SSE = FALSE, FLAG_SSI = FALSE),
                                modelling_specs = list(n_reps = 10, n_mcmc = 1000)){
  
  #1. RUN MCMC 
  #RUN_MCMC_REPS(epidemic_data, root_folder, model_type = model_type, modelling_specs = modelling_specs)
  
  #2. GET SUMMARY STATS
  GET_SUM_STATS_TOTAL(epidemic_data, root_folder, model_type = model_type)
  
  #3. GET P VALUES TOTAL
  GET_P_VALUES_TOTAL(root_folder, modelling_specs$n_reps)
  
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
                          modelling_specs = modelling_specs){
  
  'Run mcmc for n_reps iterations and save' 
  
  for(rep in 1:modelling_specs$n_reps) {
    
    #Folder 
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(root_folder, '/rep_', rep)
    print(paste0('folder_rep = ', folder_rep))
    
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
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_output, file = paste0(folder_rep, '/mcmc_output_rep_', rep, '.rds' ))
  }
}

#' Get and save the summary statistics of the real and replicated data  
#'
#'  Get and save the summary statistics of the real and replicated data. The data is replicated using a thinned sample of the MCMC output.
#'  The summary statistics of each replicated dataset are calculated and stored together in a dataframe. 
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param root_folder root folder location in which to store the MCMC results for \code{"n_reps"}
#' @param model_type A list of indicators of type boolean to indicate the type of epidemic model to fit and corresponding MCMC sampler to run;
#' \itemize{
#'   \item \code{"FLAG_BASE"} - if TRUE the baseline model MCMC sampler is implemented 
#'   \item \code{"FLAG_SSE"} - if TRUE the super-spreading events (SSE) model MCMC sampler is implemented 
#'   \item \code{"FLAG_SSI"} - if TRUE the super-spreading individuals (SSI) model MCMC sampler is implemented 
#' @param modelling_specs A list of specifications;
#' \itemize{
#'   \item \code{"n_reps"} - Number of repetitions of the MCMC sampler. Multiple runs required for model criticism
#'   \item \code{"n_mcmc"} - Number of iterations of the mcmc sampler }
#'   \item \code{"burn_in_pc"} - Proportion of mcmc samples to remove at the start as burn-in
#'   \item \code{"thinning_factor"}  - The data is replicated using a thinned sample of the MCMC output

#' @return saves the \code{"df_ss_true_rep_"} and \code{"df_summary_stats"} int the \code{"root_foler"}/\code{"rep"} location
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#'GET_SUM_STATS_TOTAL(epidemic_data, root_folder)
#' 
GET_SUM_STATS_TOTAL <- function(epidemic_data, root_folder,
                                 model_type = list(FLAG_BASE = TRUE, FLAG_SSE = FALSE, FLAG_SSI = FALSE),
                                 modelling_specs = list(n_reps = 10, n_mcmc = 10000, thinning_factor = 10, 
                                                        burn_in_size = 0.01, FLAG_THIN = TRUE)) { 
  
  'Get summary stats and p valuess for all mcmc reps'
  
  #THINNED MCMC
  if(modelling_specs$FLAG_THIN){
    mcmc_vec_size = modelling_specs$n_mcmc/modelling_specs$thinning_factor; print(paste0('thinned mcmc vec size = ', mcmc_vec_size))
  } else {
    mcmc_vec_size = modelling_specs$n_mcmc
  }
  
  #BURN-IN
  burn_in = modelling_specs$burn_in_size*modelling_specs$n_mcmc
  num_days = length(epidemic_data)
  
  #REPS FOR SUMMARY STATS
  for(rep in 1:modelling_specs$n_reps) {
    
    #REP
    print(paste0('rep = ', rep))
    folder_rep = paste0(root_folder, "/rep_", rep, '/')
    print(paste0('folder_rep = ', folder_rep))
    
    #MCMC
    mcmc_output <- readRDS(paste0(folder_rep, 'mcmc_output_rep_', rep, '.rds'))
    print('passed 1')
    
    #GET SUMMARY STATS (TRUE)
    df_ss_true = GET_SUMMARY_STATS(epidemic_data) 
    saveRDS(df_ss_true, file = paste0(folder_rep, 'df_ss_true_rep_', rep, '.rds' ))
    print('passed 2')
    
    #REPLICATED DATA (THINNED)
    for(i in seq(burn_in, mcmc_vec_size, by = modelling_specs$thinning_factor)){
      print(paste0('i = ', i))
      if (model_type$FLAG_BASE) {
        R0 = mcmc_output$r0_vec[i]
        sim_data = SIMULATE_BASELINE_EPIDEMIC(R0, num_days = num_days)
        print(sim_data)
        print('passed 3')
        
      }  else if (model_type$FLAG_SSE) {
        
        sim_data = SIMULATE_SS_EVENTS(mcmc_output$alpha_vec[i], mcmc_output$beta_vec[i], mcmc_output$gamma_vec[i])

      } else if (model_type$FLAG_SSI) {

        sim_data = SSI_MCMC_ADAPTIVE(mcmc_output$a_vec[i], mcmc_output$b_vec[i], mcmc_output$c_vec[i])
      }

      #SAVE DATA
      folder_rep_data = paste0(folder_rep, '/replicated_data/')
      ifelse(!dir.exists(file.path(folder_rep_data)), dir.create(file.path(folder_rep_data), recursive = TRUE), FALSE)
      
      saveRDS(sim_data, file = paste0(folder_rep_data, 'sim_data_iter_', i, '.rds' ))
      print('passed 4')
      
      #SUMMARY STATISTICS
      if (i == burn_in) { 
        df_summary_stats = GET_SUMMARY_STATS(sim_data)
        list_ss_iters = c(i)
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = GET_SUMMARY_STATS(sim_data, FLAG_DF_CREATE = FALSE)
        list_ss_iters = c(list_ss_iters, i)
      }
    }
    
    #SAVE SUMMARY STATISTICS + ITERATIONS
    saveRDS(df_summary_stats, file = paste0(folder_rep, '/df_replicated_summary_stats_', rep, ".rds"))
    saveRDS(list_ss_iters, file = paste0(folder_rep, '/list_ss_iters_i', rep, '.rds'))
    print('passed 5')
  }
}

#' Get and save the posterior predictive p value (ppp-value) for each summary statistic  
#'
#'  Get and save the posterior predictive p value (ppp-value) for each summary statistic across all reps. The ppp-value of each summary statistic is calculated and stored together in a dataframe. 
#'
#' @param epidemic_data data from the epidemic, namely daily infection counts
#' @param root_folder root folder location in which to store the MCMC results for \code{"n_reps"}
#' @param n_reps Total number of repetitions of the MCMC sampler ran. Ppp-values calculated on the aggregate of all reps. 
#' @return A Dataframe of ppp-values for all the summary statistics \code{"df_p_values"}. Also saves the dataframe in the \code{"root_foler"} location
#' @export
#'
#' @author Hannah Craddock, Xavier Didelot, Simon Spencer
#'
#' @examples
#'
#'GET_P_VALUES_TOTAL(root_folder, n_reps)
#' 
GET_P_VALUES_TOTAL <- function(root_folder, n_reps){
  
  for(rep in 1:n_reps) {
    
    #FOLDER REP
    print(paste0('rep = ', rep))
    folder_rep = paste0(root_folder, "/rep_", rep, '/')
    print(paste0('folder_rep = ', folder_rep))

    #GET TRUE SUMMARY STATISTICS 
    df_true_ss = readRDS(folder_rep, 'df_ss_true_rep_', rep, '.rds' )
    #df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    
    #GET REPLICATED SUMMARY STATISTICS 
    df_summary_stats_rep <- readRDS(paste0(folder_rep, '/df_summary_stats_', rep, '.rds' ))
    
    #GET P VALUES
    list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_p_vals, file = paste0(folder_rep, '/list_p_vals_rep', rep, ".rds"))
    
    #list_all_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values_list(df_summary_stats_rep[,x], df_true_ss[,x]))
    #saveRDS(list_all_p_vals, file = paste0(folder_rep, '/list_all_p_vals_rep_', rep, ".rds"))
    
    #DATAFRAME OF P VALUES 
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
                               med_dif_normIII =  list_p_vals[20])
      print(paste0('df_p_values', df_p_values))
    } else {
      df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
    }
  }
  
  #DATAFRAME
  df_p_values = as.data.frame(df_p_values)
  
  #SAVE DATAFRAME
  saveRDS(df_p_values, file = paste0(root_folder, '/total_p_values_iter_', iter, '.rds' ))
  
  #RETURN P VALUES
  df_p_values
  
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
#' @param epidemic_data Data for which the summary statistics are calculated
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
GET_SUMMARY_STATS <- function(epidemic_data, FLAG_DF_CREATE = TRUE){
  
  #Initialise
  mid_point = length(epidemic_data)/2; end = length(epidemic_data)
  
  if (FLAG_DF_CREATE){
    
    #Data-frame of summary statistics (20)
    summary_stats_results = data.frame( 
      
      sum_infects = sum(epidemic_data),
      sum_1st_half  = sum(epidemic_data[1:mid_point]),
      sum_2nd_half =  sum(epidemic_data[mid_point+1:end]),
      
      median_infect = median(epidemic_data),
      max_infect = max(epidemic_data),
      sd_infect = sd(epidemic_data),
      
      infect_q_75 = quantile(epidemic_data)[4][1][1],
      infect_q_87_5 = quantile(epidemic_data, probs = seq(0, 1, 0.125))[8][1][1],
      
      #Differences
      max_dif = max((diff(epidemic_data))), #Change from absolute difference
      med_dif = median(diff(epidemic_data)),
      max_dif2nd = max(diff(diff(epidemic_data))),
      med_dif2nd = median(diff(diff(epidemic_data))),
      
      #Norm Differences
      max_dif_normI = max(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]+1)),
      max_dif_normII = max(diff(epidemic_data)/(rollapply(epidemic_data, 2, mean, by = 1)+1)), #mean of consecutive counts 
      max_dif_normIII = max(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]/(epidemic_data[2:length(epidemic_data)]+1)+1)), #ratio of consecutive counts
      
      max_dif2nd_I = max(diff(diff(epidemic_data))/(epidemic_data[1:length(epidemic_data)-1]+1)),
      max_dif_2ndII = max(diff(diff(epidemic_data))/(rollapply(epidemic_data, 2, mean, by = 1)+1)),
      
      med_dif_normI = median(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]+1)),
      med_dif_normII = median(diff(epidemic_data)/(rollapply(epidemic_data, 2, mean, by = 1) +1)),
      med_dif_normIII = median(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]/(epidemic_data[2:length(epidemic_data)]+1)+1))
      
    )
    
  } else {
    
    #List if Dataframe already created
    summary_stats_results = list(sum(epidemic_data), sum(epidemic_data[1:mid_point]), sum(epidemic_data[mid_point +1:end]),
                                 median(epidemic_data), max(epidemic_data), sd(epidemic_data),
                                 quantile(epidemic_data)[4][1][1], quantile(epidemic_data, probs = seq(0, 1, 0.125))[8][1][1],
                                 max(diff(epidemic_data)), median(diff(epidemic_data)), max(diff(diff(epidemic_data))), median(diff(diff(epidemic_data))),
                                 max(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]+1)), max(diff(epidemic_data)/(rollapply(epidemic_data, 2, mean, by = 1)+1)),
                                 max(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]/(epidemic_data[2:length(epidemic_data)]+1)+1)),
                                 max(diff(diff(epidemic_data))/(epidemic_data[1:length(epidemic_data)-1]+1)),  max(diff(diff(epidemic_data))/(rollapply(epidemic_data, 2, mean, by = 1)+1)),
                                 median(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]+1)), median(diff(epidemic_data)/(rollapply(epidemic_data, 2, mean, by = 1)+1)),
                                 median(diff(epidemic_data)/(epidemic_data[1:length(epidemic_data)-1]/(epidemic_data[2:length(epidemic_data)]+1)+1))
                                 
    )
  }
  
  summary_stats_results
  
}