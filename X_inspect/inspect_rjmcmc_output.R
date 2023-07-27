#Inspect rjmcmc output

run = 1
OUTPUT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/rjmcmc_base1"
CURRENT_OUTPUT_FOLDER = paste0(OUTPUT_FOLDER, '/run_', run)

#RJMCMC Output
rj_out1 = readRDS(file = paste0(CURRENT_OUTPUT_FOLDER, '/rjmcmc', run, '.rds' ))
beta_vec = rj_out1$beta_vec
  
#GET BAYES FACTORS
beta_pc0 = (length(which(beta_vec == 0)))/length(beta_vec) 
bayes_factor = beta_pc0/(1-beta_pc0); bayes_factor = round(bayes_factor, 4)

plot.ts(beta_vec)