#*********************************************************
#*
#* GET POSTERIOR MODEL PROBABILITIES
#* 
#***********************************************************
#' @export 
# POSTERIOR MODEL PROB OF ONE MODEL
GET_POSTERIOR_MODEL_PROB <- function(probs_models, row_model_ev){
  
  'Given one model evidence result per model'
  
  #PARAMS
  num_models = length(row_model_ev)
  vec_model_diffs = c() 
  
  for(i in 1:(num_models-1)){
    #print(i); print(probs_models[[1]]); print(log_model_evidence[[1]])
    
    vec_model_diffs[i] = exp(log(probs_models[i+1]) + row_model_ev[i+1] -
                               log(probs_models[1]) - row_model_ev[1]) #Matches formula?
  }
  
  posterior_prob =  1/(1 + sum(vec_model_diffs))
  
  print(paste0('posterior_prob= ', posterior_prob))
  
  return(posterior_prob)
}

#VECTOR OF MODELS POSTERIOR PROBALITIES
GET_MODELS_POSTERIOR_PROBS <- function(list_log_mod_ev, num_models = 5, 
                                       probs_models = c(0.5, rep(0.5/4, 4)), FLAG_BASELINE = FALSE){
  
  'GET POSTERIOR MODEL PROBABILITY FOR EACH MODEL'
  #PARAM
  vec_post_probs = vector('numeric', length = num_models)
  
  #PRIOR PROBS: ADJUST TO ACCOUNT FOR 50% WEIGHT ON BASELINE MODEL
  if(FLAG_BASELINE){
    probs_models[1] = 0.5
  } else {
    probs_models[1] = 0.5/4
    probs_models[2] = 0.5
  }
  print(probs_models)
  
  for (i in 1:num_models){
    
    if (i == 1){
      FLAG_BASELINE = TRUE
      probs_models[1] = 0.5
    }
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(probs_models, list_log_mod_ev[i])
  }
  
  return(vec_post_probs)
  
}

#GET AGGREGATE OF POSTERIOR MODEL PROBABILITIES
GET_AGG_POSTERIOR_PROBABILITES <- function(num_models, FLAG_BASELINE = FALSE,
                                   list_log_mod_evid){ 
  
  'Get posterior probabilites for mulitple P_hat reps'
  
  #CREATE MATRIX (ROWS ARE ITERATIONS, COLUMNS ARE MODELS)
  num_reps = length(list_log_mod_evid[[1]])
  matrix_model_ev = matrix(nrow = num_reps, ncol = num_models)
  print(dim(matrix_model_ev))
  vec_post_probs = c();
  idx_sample = sample(1:num_reps, num_reps)
  
  #PROBABILITY
  probs_models = c()
  for (i in 1:num_models){
    probs_models = c(probs_models, 0.5/num_models)
    matrix_model_ev[,i] = list_log_mod_evid[[i]]
  }
  
  #ADJUST PROBS_MODELS TO ACCOUNT FOR 50% WEIGHT ON BASELINE MODEL
  if(FLAG_BASELINE){
    probs_models[1] = 0.5
  } else {
    probs_models[2] = 0.5
  }
  print(probs_models)
  
  for (i in 1:num_reps){
    print(paste0('rep = ', i))
    j = idx_sample[i]; print(j)
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(probs_models, matrix_model_ev[j,])
  }
  
  return(vec_post_probs)
}
