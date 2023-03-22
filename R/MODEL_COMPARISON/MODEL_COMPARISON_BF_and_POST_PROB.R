#MODEL COMPARISON

#*********************************************************
#* GET BAYES FACTORS
#*******************

GET_BAYES_FACTORS <- function(loglike_vec1, loglike_vec2){
  
  'Get Bayes factor via ratio of the model evidence'
  
  #bayes_factor = MODEL_EVIDENCE(loglike_vec1)/MODEL_EVIDENCE(loglike_vec2) 
  log_bf = LOG_MODEL_EVIDENCE(loglike_vec1)  - LOG_MODEL_EVIDENCE(loglike_vec2)
  bayes_factor = exp(log_bf)
  
  return(bayes_factor)
}


#*********************
#* 4. GET BAYES FACTORS
#*********************
GET_LOG_BAYES_FACTORS <- function (list_log_phat_mod1, list_log_phat_mod2){
  
  log_bf = list_log_phat_mod1 - list_log_phat_mod2
  
  return(log_bf)
  
}

#*********************************************************
#*
#* 2. GET POSTERIOR MODEL PROBABILITIES
#* 
#***********************************************************
GET_POSTERIOR_MODEL_PROB <- function(num_models = 3, 
                                     probs_models = list(prob_mech1 = 0.25, prob_mech2 = 0.5, prob_mech3 = 0.25), #mech = mechanism; baseline (0.5) or sse (0.25)
                                     log_model_evidence = list(mod1 = mod1,
                                                      mod2 = mod2, mod3 = mod3)){ 
  
  #PARAMS
  vec_model_diffs = c()
  
  for(i in 1:num_models-1){
    vec_model_diffs[i] = exp(log(probs_models[[i+1]]) + log_model_evidence[[i+1]] -
                               log(probs_models[[1]]) - log_model_evidence[[1]])
  }
  
  print(paste0('sum(vec_model_diffs) = ', sum(vec_model_diffs)))
  
  posterior_prob =  1/(1 + sum(vec_model_diffs))
  
  print(paste0('posterior_prob= ', posterior_prob))
  
  return(posterior_prob)
}

#GET_AGGREGATE_POSTERIOR_MODEL_PROB
GET_AGG_POSTERIOR_PROB <- function(num_models = 3, 
                                   probs_models = list(prob_mech1 = 0.25, prob_mech2 = 0.5, prob_mech3 = 0.25),
                                   FLAG_EQUAL_PROB = TRUE,
                                   list_log_mod_evid = list(mod1 = mod1,
                                                         mod2 = mod2, mod3 = mod3)){ 
  'Get posterior probabilites for mulitple P_hat reps'
  
  #FOR EACH REP
  num_reps = length(list_log_mod_evid$mod1)
  vec_post_probs = c()
  #PROBABILITY
  if(FLAG_EQUAL_PROB){
    probs_models = list(prob_mech1 = 1/num_models, prob_mech2 = 1/num_models, prob_mech3 = 1/num_models) 
    print(probs_models)
  }
  
  for (i in 1:num_reps){
    print(paste0('rep = ', i))
    vec_post_probs[i] = GET_POSTERIOR_MODEL_PROB(num_models = 3, 
                                                 probs_models = probs_models,
                                                 log_model_evidence = list(mod1 = list_log_mod_evid$mod1[i],
                                                                  mod2 = list_log_mod_evid$mod2[i],
                                                                  mod3 = list_log_mod_evid$mod3[i]))
    
  }
  
  return(vec_post_probs)
  
}