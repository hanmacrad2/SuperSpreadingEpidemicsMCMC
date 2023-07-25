#INSPECT SD OF MODEL EVIDENCE SD OF SSIB MODEL

REPEAT_MODEL_EV <- function(reps = 5){
  
  vec_sd = vector(length = reps)
  for (i in 1:reps){
    
    print('*********************')
    print(paste0('i = ', i))
    model_ev_ssib10 = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                   FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                       SSIB = TRUE, SSIR = FALSE))
    print(mean(model_ev_ssib10))
    print(paste0('sd(model_ev_ssib10)',  sd(model_ev_ssib10)))
    vec_sd[i] = sd(model_ev_ssib10)
  } 
  
  return(vec_sd)
}

mod_ev_orig = c(2.22, 4.11, 4.03, 4.86, 5.71)
mean(mod_ev_orig)

model_ev_ac = c(3.7, 4.25, 5.77, 5.83, 3.35)
mean(model_ev_ac)

#Rep 4, beta = 10
sd_10 = REPEAT_MODEL_EV()

#Rep 6, beta = 0.1
sd_10 = REPEAT_MODEL_EV()

#Rep 7, num is = 10000
sd_10 = REPEAT_MODEL_EV()

sd_ac_all = REPEAT_MODEL_EV()
mean(sd_ac_all)

sd_ac_b01 = REPEAT_MODEL_EV()
mean(sd_ac_b01)

sd_ac_b05 = REPEAT_MODEL_EV()
mean(sd_ac_b05)

#*****************
#IS SAMPS = 10K
sd_ac_b01_10k = REPEAT_MODEL_EV()
mean(sd_ac_b01_10k)

sd_ac_b005_10k = REPEAT_MODEL_EV()
mean(sd_ac_b005_10k)

sd_ac_b01_10k_no_ac = REPEAT_MODEL_EV()
mean(sd_ac_b01_10k_no_ac)

sd_ac_b005_10k_no_ac = REPEAT_MODEL_EV()
  mean(sd_ac_b005_10k_no_ac)
  