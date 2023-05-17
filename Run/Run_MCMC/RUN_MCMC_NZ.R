#NEW MCMC SCRIPT MULTI TIMES SCRIPT_3

#RUN MULTIPLE MCMC ITERATIONS
library(SuperSpreadingEpidemicsMCMC)
library(MASS)
ls("package:SuperSpreadingEpidemicsMCMC")

#FOLDER RESULTS
data_nz = 'NZ_DATA_WAIT_21_SUBSET_I/'
data_nz = 'NZ_DATA_WAIT_21/'
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, data_nz)
print(OUTER_FOLDER)
create_folder(OUTER_FOLDER)
run_number = 1

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'CM_08_21_SUB_1/')
print(OUTER_FOLDER)
create_folder(OUTER_FOLDER)
run_number = 1

#***********************
# EPIDEMIC DATA -- NEW ZEALAND
#**********************
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
filename =  'data_waitemata_aug_21.csv'

data_file_wait_21 = read.csv(paste0(DATA_FOLDER, filename))
data_wait_08_21 = data_file_wait_21$Cases
plot.ts(data_wait_08_21, ylab = 'Infection count', main = 'Waitemata NZ, August 2021')
saveRDS(data_wait_08_21, file = paste0(OUTER_FOLDER, 'data_wait_08_21.rds'))

#SUBSET
data_wait_08_21_sub1 = data_wait_08_21[2:11]
plot.ts(data_wait_08_21_sub1, ylab = 'Infection count', main = 'Waitemata NZ, August 2021')

#REMOVE LEADING ZEROES FOR MCMC
#data_wait_08_21 = c(0, 0, data_wait_08_21)
data_wait_08_21 = data_wait_08_21[1:length(data_wait_08_21)]
data_wait_08_21 = data_wait_08_21[2:length(data_wait_08_21)]

#********************
#DATA: COUNTIES MANUKAAU
#*********************
filename =  'data_cm_tot_aug_21.csv'

data_file_cm_08_21 = read.csv(paste0(DATA_FOLDER, filename))
data_wait_cm_08_21 = data_file_cm_08_21$Cases
plot.ts(data_wait_cm_08_21, ylab = 'Infection count', main = 'COUNTIES MANAUKA NZ, August 2021')


data_file_cm_08_21_sub1 = data_wait_cm_08_21[5:21]
plot.ts(data_file_cm_08_21_sub1)
data_file_cm_08_21_sub2 = data_wait_cm_08_21[25:47]
plot.ts(data_file_cm_08_21_sub2)

#***********************
# 2. RUN BASELINE MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_file_cm_08_21_sub1, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))


#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_file_cm_08_21_sub1, OUTER_FOLDER, run_number = 2, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_file_cm_08_21_sub1, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN SSIR MCMC
#**********************
run = 4
RUN_MCMC_MULTIPLE_TIMES(data_wait_08_21_sub1, OUTER_FOLDER, run_number = run, n_repeats = 10, n_mcmc = 50000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIR = TRUE))

#***********
#DATA 2
#***********
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'CM_08_21_SUB_2/')
print(OUTER_FOLDER)

#BASE
RUN_MCMC_MULTIPLE_TIMES(data_file_cm_08_21_sub2, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))


#***********************
# 3. RUN SSNB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_file_cm_08_21_sub2, OUTER_FOLDER, run_number = 2, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE, SSNB = TRUE,
                                            SSIB = FALSE, SSIC = FALSE))

#***********************
# 2. RUN SSEB MCMC
#**********************
RUN_MCMC_MULTIPLE_TIMES(data_file_cm_08_21_sub2, OUTER_FOLDER, run_number = run_number, n_repeats = 50, n_mcmc = 30000,
                        FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE, SSNB = FALSE,
                                            SSIB = FALSE, SSIC = FALSE))



#**************
#INSPECT MCMC OUTPUT
#**************

i = 1
model_type = 'SSIR'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/')
mcmc_ssir =  readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))

i = 1
model_type = 'SSEB'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/')
mcmc_sseb =  readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))

#PLOT
PLOT_SSB_MCMC_REAL_DATA(data_wait_08_21_sub1, mcmc_sseb, 30000)

#INSPECT
mcmc_output = MCMC_INFER_BASELINE(data_wait_08_21, n_mcmc = 10) #00)

mcmc_output = MCMC_INFER_SSEB(data_wait_08_21, n_mcmc = 1000)   


for (t in 2:num_days) {
  
  #ETA (t-1)
  eta_vec[t-1] <- rgamma(1, shape = x[t-1]*k, scale = R0X/k) #Draw eta from previous time step
  #INFECTIVITY
  infectivity = rev(prob_infect[1:t]) 
  #POISSON; OFFSPRINT DISTRIBUTION
  total_rate = sum(eta_vec*infectivity) #DOT PRODUCT
  x[t] = rpois(1, total_rate)
  
}
