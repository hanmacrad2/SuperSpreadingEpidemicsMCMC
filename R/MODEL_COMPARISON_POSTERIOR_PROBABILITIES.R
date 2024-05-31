#*********************************************************
#*
#* GET POSTERIOR MODEL PROBABILITIES
#* 
#***********************************************************
#' @export 
GET_MODEL_EVIDENCE <- function(list_mcmc, epidemic_data){
  
  vec_mod_ev = c()
  
  #***************************
  #1. BASELINE MODEL
  print('BASELINE:')
  mcmc_baseline = list_mcmc$mcmc_baseline #MCMC_INFER_BASELINE(epidemic_data, n_mcmc)
  #MODEL EVIDENCE 
  mod_ev_base = GET_LOG_MODEL_EVIDENCE_BASELINE(mcmc_baseline$r0_vec, epidemic_data)
  print(mod_ev_base)
  vec_mod_ev = c(vec_mod_ev, mod_ev_base)
  
  #SSE 
  print('SSE:')
  mcmc_sse = list_mcmc$mcmc_sse
  #MODEL EVIDENCE
  mcmc_samples =  mcmc_sse$sse_params_matrix 
  FLAGS_MODELS = GET_FLAGS_MODELS(SSE = TRUE) 
  mod_ev_sse = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data, FLAGS_MODELS = FLAGS_MODELS) 
  print(mod_ev_sse)
  vec_mod_ev = c(vec_mod_ev, mod_ev_sse)
  
  #SSI
  print('SSI:')
  mcmc_ssi = list_mcmc$mcmc_ssi
  #MODEL EVIDENCE
  mod_ev_ssi = GET_LOG_MODEL_EVIDENCE_SSI(mcmc_ssi, epidemic_data) 
  print(mod_ev_ssi)
  vec_mod_ev = c(vec_mod_ev, mod_ev_ssi)
  
  #SSEB
  print('SSEB:')
  mcmc_sseb = list_mcmc$mcmc_sseb
  #MODEL EVIDENCE
  mcmc_samples =  matrix(c(mcmc_sseb$alpha_vec, mcmc_sseb$r0_vec,
                           mcmc_sseb$beta_vec), ncol = 3)
  FLAGS_MODELS = GET_FLAGS_MODELS(SSEB = TRUE) 
  mod_ev_sseb = GET_LOG_MODEL_EVIDENCE(mcmc_samples, epidemic_data,
                                       FLAGS_MODELS = FLAGS_MODELS) 
  print(mod_ev_sseb)
  vec_mod_ev = c(vec_mod_ev, mod_ev_sseb)
  
  #SSIB
  print('SSIB:')
  mcmc_ssib = list_mcmc$mcmc_ssib
  #MODEL EVIDENCE
  mod_ev_ssib = GET_LOG_MODEL_EVIDENCE_SSIB(mcmc_ssib, epidemic_data)
  print(mod_ev_ssib)
  vec_mod_ev = c(vec_mod_ev, mod_ev_ssib)
  
  print('vec mod ev'); print(vec_mod_ev)
  
  return(vec_mod_ev)
}




GET_MODELS_POSTERIOR_PROBS <- function(vec_log_mod_ev, num_models = 5, 
                                       prior_probs_models = c(0.5, rep(0.5/4, 4))){
  
  'Get posterior probabilities for each model vs other 4 models in vec_log_mod_ev'
  # PARAM
  post_probs_vec = vector('numeric', length = num_models)
  
  for (i in 1:length(vec_log_mod_ev)){
    
    #PARAMS
    mod_ev_i = vec_log_mod_ev[i]; prior_mod_i = prior_probs_models[i]
    vec_model_diffs = vector('numeric', length = num_models - 1)
    #print(paste0('i: ', i, 'mod_ev_i:', mod_ev_i))
    
    #Loop over other models
    for (j in seq_along(vec_log_mod_ev)) {
      if (j != i) {
        #browser()
        vec_model_diffs[j] = exp(vec_log_mod_ev[j] + log(prior_probs_models[j])-
                                   mod_ev_i - log(prior_mod_i))
      }
      
    }
    #Posterior_prob
    post_probs_vec[i] =  1/(1 + sum(vec_model_diffs))
  }
  
  print(post_probs_vec)
  return(post_probs_vec)
}


GET_MODEL_COMP_DF_PLOT <- function(df_mod_ev_plot){
  
  'Return: df_post_probs'
  
  df_mod_ev_plot <- df_mod_ev_plot[, !colnames(df_mod_ev_plot) %in% "tot_infs"]
  df_mod_ev_plot <- df_mod_ev_plot[, !colnames(df_mod_ev_plot) %in% "run_number"]
  df_mod_ev_plot <- df_mod_ev_plot[, !colnames(df_mod_ev_plot) %in% "data_sim"]
  df_mod_ev_plot <- df_mod_ev_plot[, !colnames(df_mod_ev_plot) %in% "sim_model"]
  
  
  df_post_probs = apply(df_mod_ev_plot, MARGIN=1, FUN=function(x2) GET_MODELS_POSTERIOR_PROBS(x2, prior_probs_models = rep(0.2, 5)) )
  df_post_probs = t(df_post_probs)
  
  return(df_post_probs)
  
}

#************************************************************************************************
#* GET MODEL COMPARISON POST PROBS ALL MODELS
GET_POST_PROB_HEAT_MAP_VALUES <- function(df_post_probs, model_num,
                                          FLAGS_MODELS, rounding_number = 3){
  
  model = names(FLAGS_MODELS)[model_num]
  print(paste0(model, ': SIMULATION MODEL'))
  
  #BASELINE FIRST
  get_mean_cis(df_post_probs[,1], rounding_number = rounding_number) #run_number = column 1
  mean_post_prob = get_mean(df_post_probs[,1], rounding_number = rounding_number)
  vec_means = c(mean_post_prob)
  list_post_probs = list(Baseline = mean_post_prob)
  
  for (i in 2:length(FLAGS_MODELS)){
    
    model_fitted = names(FLAGS_MODELS)[i]
    get_mean_cis(df_post_probs[,i], rounding_number = rounding_number)
    mean_post_prob = get_mean(df_post_probs[,i], rounding_number = rounding_number)
    print(paste0(model_fitted, ': FITTED MODEL, mean post prob: ', mean_post_prob))
    print(paste0('model (should be same): ', colnames(df_post_probs)[i]))
    list_post_probs[model_fitted] = mean_post_prob
    vec_means = c(vec_means, mean_post_prob)
    
  }
  
  print(list_post_probs)
  print(vec_means)
  
  list_results = list(list_post_probs = list_post_probs, vec_means = vec_means)
  
  return(list_results)
}


#****************
#* MAP : MAXIMUM A POSTERIOR

GET_COUNT_MAP_ALL_MODELS <- function(matrix_probs, model_num, FLAG_MODEL){
  
  #MODEL
  FLAGS_MODELS = list(Baseline = FALSE, SSE = TRUE, SSI = FALSE,
                      SSEB = FALSE, SSIB = FALSE) 
  model = names(FLAG_MODEL)[which(unlist(FLAG_MODEL))]
  model_index = which(unlist(FLAGS_MODELS))[[1]]
  print(model)
  
  for (i in 1:length(FLAGS_MODEL)){
    is_max_column <- function(row, model_num) {
      if (row[model_num] == max(row)) {
        return(1)
      } else {
        return(0)
      }
    }
    
    counts_model_map <- apply(matrix_probs, 1, function(row) is_max_column(row, model_num))
  }
  
  # Function to determine if column 1 is the maximum value in the row

  
  # Print the sum of counts
  print(paste('count_model:', sum(counts_model_map)))
}


GET_COUNT_MAP <- function (matrix_probs, model_num){
  
  # Function to determine if column 1 is the maximum value in the row
  is_max_column <- function(row, model_num) {
    if (row[model_num] == max(row)) {
      return(1)
    } else {
      return(0)
    }
  }
  
  counts_model_map <- apply(matrix_probs, 1, function(row) is_max_column(row, model_num))
  
  # Print the sum of counts
  print(paste('count_model:', sum(counts_model_map)))
}
