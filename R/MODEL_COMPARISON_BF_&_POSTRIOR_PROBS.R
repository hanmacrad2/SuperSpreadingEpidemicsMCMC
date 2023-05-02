#MODEL COMPARISON

#*********************************************************
#* GET BAYES FACTORS
#*******************

# GET_BAYES_FACTORS <- function(loglike_vec1, loglike_vec2){
#   
#   'Get Bayes factor via ratio of the model evidence'
#   
#   #bayes_factor = MODEL_EVIDENCE(loglike_vec1)/MODEL_EVIDENCE(loglike_vec2) 
#   log_bf = LOG_MODEL_EVIDENCE(loglike_vec1)  - LOG_MODEL_EVIDENCE(loglike_vec2)
#   bayes_factor = exp(log_bf)
#   
#   return(bayes_factor)
# }
#' @export
GET_BAYES_FACTORS <- function(list_log_phat_mod1, list_log_phat_mod2){
  
  log_bf = list_log_phat_mod1 - list_log_phat_mod2
  bayes_factor = exp(log_bf)
  
  return(bayes_factor)
  
}

#*********************
#* 4. GET BAYES FACTORS
#*********************
#' @export
GET_LOG_BAYES_FACTORS <- function (list_log_phat_mod1, list_log_phat_mod2){
  
  log_bf = list_log_phat_mod1 - list_log_phat_mod2
  
  return(log_bf)
  
}

#*********************************************************
#*
#* 2. GET POSTERIOR MODEL PROBABILITIES
#* 
#***********************************************************
#' @export
GET_POSTERIOR_MODEL_PROB <- function(log_model_evidence, probs_models){
  
  #PARAMS
  num_models = length(log_model_evidence)
  print(paste0('num models = ', num_models))
  vec_model_diffs = c() 
  
  for(i in 1:(num_models-1)){
    #print(i); print(probs_models[[1]]); print(log_model_evidence[[1]])
    
    vec_model_diffs[i] = exp(log(probs_models[[i+1]]) + log_model_evidence[[i+1]] -
                               log(probs_models[[1]]) - log_model_evidence[[1]]) #Matches formula?
  }
  
  print(paste0('sum(vec_model_diffs) = ', sum(vec_model_diffs)))
  
  posterior_prob =  1/(1 + sum(vec_model_diffs))
  
  print(paste0('posterior_prob= ', posterior_prob))
  
  return(posterior_prob)
}

#GET_AGGREGATE_POSTERIOR_MODEL_PROB
#' @export
GET_AGG_POSTERIOR_PROB <- function(num_models = 3, FLAG_EQUAL_PROB = FALSE, probs_models =
                                     list(prob1 = 0.5, prob2 = 0.25, prob3 = 0.25), list_log_mod_evid){ 
  'Get posterior probabilites for mulitple P_hat reps'
  
  #FOR EACH REP
  num_reps = length(list_log_mod_evid$mod1)
  vec_post_probs = c();
  idx_sample = sample(1:num_reps, num_reps)
  
  #PROBABILITY
  if(FLAG_EQUAL_PROB){
    probs_models = list(prob_mech1 = 1/num_models, prob_mech2 = 1/num_models, prob_mech3 = 1/num_models) 
  }
  print(probs_models)
  
  for (i in 1:num_reps){
    print(paste0('rep = ', i))
    j = idx_sample[i]; print(j)
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(log_model_evidence = list(mod1 = list_log_mod_evid$mod1[j],
                                                                  mod2 = list_log_mod_evid$mod2[j],
                                                                  mod3 = list_log_mod_evid$mod3[j]), 
                                                 probs_models)
    
  }
  
  return(vec_post_probs)
  
}

#MODELS
#' @export
GET_AGG_POSTERIOR_PROB_II <- function(num_models = 2, FLAG_EQUAL_PROB = TRUE,                                               #probs_models = list(prob_mech1 = 0.25, prob_mech2 = 0.5, prob_mech3 = 0.25)
                                   list_log_mod_evid){ 
  'Get posterior probabilites for mulitple P_hat reps'
  
  #FOR EACH REP
  num_reps = length(list_log_mod_evid$mod1)
  vec_post_probs = c();
  idx_sample = sample(1:num_reps, num_reps)
  
  #PROBABILITY
  if(FLAG_EQUAL_PROB){
    probs_models = list(prob1 = 1/num_models, prob2 = 1/num_models) 
  }
  print(probs_models)
  
  for (i in 1:num_reps){
    #print(paste0('rep = ', i))
    j = idx_sample[i]; #print(j)
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(log_model_evidence = list(mod1 = list_log_mod_evid$mod1[j],
                                                                           mod2 = list_log_mod_evid$mod2[j]),
                                                 probs_models)
    
  }
  
  return(vec_post_probs)
  
}