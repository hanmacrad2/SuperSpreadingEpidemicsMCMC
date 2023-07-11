#**********************************************
#GET MODEL EVIDENCE VIA IMPORTANCE SAMPLING (2018 Paper)
#***********************
library(SuperSpreadingEpidemicsMCMC)
library(mvtnorm)

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/BASELINE_DATA/DATA_BASELINE_1/"

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
#OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/PART_2/SSEB_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21_SUBSET_I/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21_SUBSET_II/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21/"

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/BASELINE_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/BASELINE_DATA/DATA_BASELINE_2/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/BASELINE_DATA/MISSING/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/SSIB_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/MOCK_DATA/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/MOCK_DATA/MOCK_DATA_4_DAYS/"

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/SSIB_DATA/DATA_SSIB_2/"
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/SSIB_DATA/DATA_SSIB_4/"

run = 1 
run = 2
n_repeats = 10

#***********************
# 1. DATA 
#**********************
EPI_DATA = data_baseline
EPI_DATA = data_ssib4 #data_ssib3
file_name = 'data_ssib3.rds'
saveRDS(data_ssib3, paste0(OUTER_FOLDER, file_name)) 

#BASE_DATA = FALSE; SSEB_DATA = TRUE; NZ_DATA = FALSE
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"

#BASE DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/BASE_DATA/"
file_name = 'epi_data_base_1.rds'
data_baseline = readRDS(file = paste0(OUTER_FOLDER, file_name))
plot.ts(data_baseline)

file_name = 'data_baseline2.rds'
data_baseline2 = readRDS(file = paste0(DATA_FOLDER, file_name))
plot.ts(data_baseline2)

#SSEB DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/SSEB_DATA/"
file_name = "epi_data_sseb_1.rds"
data_sseb = readRDS(file = paste0(OUTER_FOLDER, file_name))
plot.ts(data_sseb)

#NZ DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21/"
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
data_file_wait_21 = read.csv(paste0(DATA_FOLDER, 'data_waitemata_aug_21.csv'))
data_wait_08_21 = data_file_wait_21$Cases
plot.ts(data_wait_08_21)

#MOCK DATA
EPI_DATA = MOCK_DATA_3_DAYS

#*************************************************
#*
#* MODEL EVIDENCE ESTIMATES
#* 
#****************************************************************************
run = 1
n_repeats = 5 #10

#*************************
#1. BASELINE
#*************************
model_ev_base5 = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = FALSE))
mean(model_ev_base5)
sd(model_ev_base5) 
#PLOT_MODEL_EV_RESULTS(model_ev_base, model_type = '. Baseline w/ exp(1) prior on R0')

#*************************
#2. SSNB
#*************************
model_ev_ssnb4 = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                        SSIB = FALSE, SSIR = FALSE))

mean(model_ev_ssnb4)
sd(model_ev_ssnb4)

#*************************
#3. SSEB
#*************************
model_ev_sseb5 = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = FALSE))
#PLOT_MODEL_EV_RESULTS(model_ev_sseb)
mean(model_ev_sseb5)
sd(model_ev_sseb5)

#*************************
#5. SSIB
#*************************
model_ev_ssib4 = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                             FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                 SSIB = TRUE, SSIR = FALSE))
mean(model_ev_ssib4)
sd(model_ev_ssib4)

#*************************
#4. SSIR
#*************************
model_ev_ssir6 = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                    FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                        SSIB = FALSE, SSIR = TRUE))
#PLOT
#PLOT_MODEL_EV_RESULTS(model_ev_ssir)
mean(model_ev_ssir5)
sd(model_ev_ssir5)

#PLOT
PLOT_MODEL_EV_RESULTS(model_ev_ssib, model_type = 'SSIB',
                      data_type = 'SSIB')

#*************************
#PLOT ALL MODEL EVIDENCE
data_type = 'Mock Data'
BOX_PLOT_MODEL_EV(list_vec_results = list(BASE = model_ev_base4, SSEB = model_ev_sseb4,
                                          SSNB = model_ev_ssnb4, SSIB = model_ev_ssib4), data_type = data_type)

#***************************************************************************
#*
#* GAMMA PRIORS
#* 
#**************************************************************************

#BASELINE + GAMMA PRIOR 
run = '1_ga'
model_ev_base_gp = LOAD_MCMC_GET_MODEL_EVIDENCE(EPI_DATA, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                FLAGS_MODELS = list(BASE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                                                    SSIB = FALSE, SSIR = FALSE))
mean(model_ev_base_gp)# , na.rm = TRUE)
sd(model_ev_base_gp) #, na.rm = TRUE)
PLOT_MODEL_EV_RESULTS(model_ev_base_gp, model_type = '. Baseline w/ gamma(2, mean = R0) prior on R0')

#*************************
#SSIB + GAMMA PRIOR 
#*************************
run = '2_ga_prior'
model_ev_ssib2 = LOAD_MCMC_GET_MODEL_EVIDENCE(data_ssib2, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                             FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                                                 SSIB = TRUE, SSIR = FALSE))
mean(model_ev_ssib2)
sd(model_ev_ssib2)
PLOT_MODEL_EV_RESULTS(model_ev_ssib2, model_type = 'SSIB; prior(a) ~ gamma(mean ~ a)',
                      data_type = 'SSIB')

#*************************
#SSNB + GAMMA PRIOR ON K
#*************************
run = 'gamma_prior_k'; 
model_ev_ssnb_gak = LOAD_MCMC_GET_MODEL_EVIDENCE(data_wait_08_21_sub2, OUTER_FOLDER, run = run, n_repeats = n_repeats,
                                                       FLAGS_MODELS = list(BASE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                                                           SSIB = FALSE, SSIR = FALSE))

PLOT_MODEL_EV_RESULTS(model_ev_ssnb_gak)

#PLOT
data_type = 'Waitemata days 10-20'
par(mfrow = c(1,1))
BOX_PLOT_MODEL_EV(list_vec_results = list(ssnb_k_exp_prior = model_ev_ssnb,
                                          ssnb_k_ga_prior = model_ev_ssnb_gak),
                              title = 'SSNB Model k priors compared exp(1) vs Ga(0.001, rt = 0.001). ',
                              data_type = data_type, model = '') 

#*****************************
#*
#LOAD MODEL EVIDENCE ESTIMATES
#*
#*****************************
LOAD_MODEL_EVIDENCE <- function(model_type, run, OUTER_FOLDER){
  

  print(model_type)
  CURRENT_FOLDER = paste0(OUTER_FOLDER, toupper(model_type), '/run_', run, '/')
  list_model_ev = readRDS(file = paste0(CURRENT_FOLDER, '/model_evidence_', model_type, '_', run, '.rds'))
  print(list_model_ev)
  return(list_model_ev)
}

#LOAD 
run = 2
model_ev_base2 = LOAD_MODEL_EVIDENCE('baseline', run, OUTER_FOLDER)
model_ev_sseb2 = LOAD_MODEL_EVIDENCE('sseb', run, OUTER_FOLDER)
model_ev_ssnb2 = LOAD_MODEL_EVIDENCE('ssnb', run, OUTER_FOLDER)
model_ev_ssir2 = LOAD_MODEL_EVIDENCE('ssir', run, OUTER_FOLDER)

#SD OF RESULTS
sd(model_ev_ssir2, na.rm = TRUE) 
