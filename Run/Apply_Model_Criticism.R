#APPLY MODEL CRITICISM
library(SuperSpreadingEpidemicsMCMC)
#library(coda)

#FOLDER
ROOT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism_II/"
ITER = 'iter_1/'
MODEL_TYPE = 'baseline/'
RESULTS_FOLDER =  paste0(ROOT_FOLDER, ITER, MODEL_TYPE)
print(RESULTS_FOLDER)

#MCMC
modelling_specs = list(n_reps = 10, n_mcmc = 1000)

#APPLY
RUN_MODEL_CRITICISM(canadaX, RESULTS_FOLDER, modelling_specs = modelling_specs)

