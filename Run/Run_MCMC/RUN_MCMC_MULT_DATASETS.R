#DATASETS BASELINE
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/PART_2/SSEB_DATA/"
run = 1
OUTER_FOLDER = paste0(OUTER_FOLDER, 'run_', run, '/')
OUTER_FOLDER
create_folder(OUTER_FOLDER)

seedX = 1
set.seed(seedX)

#PARAMS
num_days = 50
n_datasets = 30

#*****************
#* 1. BASE DATA
#******************

#SETUP MATRIX DATA
R0X = 1.6; 
matrix_data_base = matrix(NA, n_datasets, num_days) 

#Simulate data
data_baseX = SIMULATE_BASELINE_EPIDEMIC(R0X)
plot.ts(data_baseX)

#Store
matrix_data_base[seedX, ] = data_baseX
matrix_data_base

seedX = seedX + 1

#*****************
#* 2. SSEB DATA
#******************
model_type = 'SSEB'
matrix_data_sseb = matrix(NA, n_datasets, num_days) 

#Simulate data
data_sseb = SIMULATE_EPI_SSEB()
plot.ts(data_sseb)

#Store
matrix_data_sseb[seedX, ] = data_sseb
matrix_data_sseb[seedX, ]

seedX = seedX + 1
set.seed(seedX)

#SAVE MATRIX
saveRDS(matrix_data_sseb, file = paste0(OUTER_FOLDER, 'matrix_datasets_', tolower(model_type), '.rds'))

#***************
#* RUN MCMC 
#***************
RUN_MCMC_MULTIPLE_DATASETS(matrix_data_sseb, OUTER_FOLDER)

#matrix_data_sseb[29,]


#Inspect results
model_type = 'SSEB'
list_is_log_ev_sseb = readRDS(paste0(OUTER_FOLDER, model_type, '/phat_ests_sseb_', run, '.rds'))
