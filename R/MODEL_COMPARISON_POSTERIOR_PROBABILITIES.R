#*********************************************************
#*
#* GET POSTERIOR MODEL PROBABILITIES
#* 
#***********************************************************
#' @export 

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


#DATAFRAME OF POSTERIOR PROBS FOR MULTIPLE SIMULATIONS/MODEL EVIDENCE RUNS
#df_post_probs = apply(df_mod_ev_plot, MARGIN=1, FUN=function(x2) GET_MODELS_POSTERIOR_PROBS(x2) )

