#****************************************************************
# DATA
#****************************************************************
#library(SuperSpreadingEpidemicsMCMC)
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/X_data/DATA_SETS/"

#BASELINE
file_name = 'data_baseline_1_5.rds'
data_baseline2 = SIMULATE_EPI_BASELINE(1.5)
plot.ts(data_baseline2, main = 'Baseline Data II, R0 = 1.5')
saveRDS(data_baseline2, paste0(DATA_FOLDER, file_name))
data_baseline = readRDS(paste0(DATA_FOLDER, file_name))

file_name = 'data_baseline2.rds'
data_baseline = readRDS(paste0(DATA_FOLDER, file_name))
plot.ts(data_baseline2)
data_base = readRDS(paste0(DATA_FOLDER, file_name))

file_name = 'data_baseline30.rds'
data_baseline30 = SIMULATE_EPI_BASELINE(2.1, num_days = 30)
data_baseline30 = readRDS(paste0(DATA_FOLDER, file_name))
EPI_DATA = data_baseline30
plot.ts(data_baseline30)
saveRDS(data_baseline30, paste0(DATA_FOLDER, file_name))

#SSNB
file_name = 'data_ssnb_1_6_2.rds'
data_ssnb2 = SIMULATE_EPI_SSNB()
plot.ts(EPI_DATA, main = 'SSE Data, R0 = 1.6')
saveRDS(data_ssnb, paste0(DATA_FOLDER, file_name))
data_ssnb = readRDS(paste0(DATA_FOLDER, file_name))

#SSIR
file_name = 'data_ssir2.rds'
data_ssir2 = SIMULATE_EPI_SSIR()
plot.ts(data_ssir$epidemic_data)
saveRDS(data_ssir, paste0(DATA_FOLDER, file_name))

data_ssir2 = readRDS(paste0(DATA_FOLDER, file_name))
data_ssir2 = data_ssir2$epidemic_data

#SSIB
file_name = 'data_ssib2.rds'
data_ssib3 = SIMULATE_EPI_SSIB(aX = 0.8)
plot.ts(data_ssib3,  main = 'SSIB DATA, R0 = 1.8', ylab = 'Daily infection count')
saveRDS(data_ssib3, paste0(DATA_FOLDER, file_name))

data_ssib2 = readRDS(paste0(DATA_FOLDER, file_name))
plot.ts(data_ssib2)
par(mfrow = c(1,1))

file_name = 'data_ssib4.rds'
data_ssib4 = SIMULATE_EPI_SSIB(num_days = 100, aX = 0.8, bX = 0.05)
plot.ts(data_ssib4, main = 'SSIB DATA, R0 = 1.3. a: 0.8, b: 0.05, c: 10', ylab = 'Daily infection count')
saveRDS(data_ssib4, paste0(DATA_FOLDER, file_name))
data_ssib4 = readRDS(paste0(DATA_FOLDER, file_name))

#MOCK DATA
EPI_DATA = MOCK_DATA_2_DAYS #= c(1,0)
 #_7_DAYS
MOCK_DATA = c(1,2,1,2,2)
MOCK_DATA = c(1,2)
MOCK_DATA_2_DAYS = c(1,0)
MOCK_DATA_3_DAYS = c(1,0,0) #c(1,2,1)
MOCK_DATA_3_DAYS = c(1,1,1)
MOCK_DATA_4_DAYS = c(1,1,0,1)
MOCK_DATA_6_DAYS = c(1,1,0,1,0,1)

MOCK_DATA_2_DAYS = c(1,0); EPI_DATA = MOCK_DATA_2_DAYS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/MOCK_DATA/MOCK_DATA_2_DAYS/"

MOCK_DATA_3_DAYS = c(1,0,0); EPI_DATA = MOCK_DATA_3_DAYS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/MOCK_DATA/MOCK_DATA_3_DAYS/"

MOCK_DATA_7_DAYS = c(1,0,0,0,1,0,1); EPI_DATA = MOCK_DATA_7_DAYS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/MOCK_DATA/MOCK_DATA_7_DAYS/"

MOCK_DATA_8_DAYS = c(1,1,1,1,1,1,1,1)
MOCK_DATA_10_DAYS = c(1,0,0,0,0,1,0,0,1,1); EPI_DATA = MOCK_DATA_10_DAYS
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/models/MOCK_DATA/MOCK_DATA_10_DAYS/"

#PLOT
plot.ts(EPI_DATA, main = 'MOCK DATA', ylab = 'Daily infection count', )

#SAVE
file_name = 'mock_data_6_days.rds'
saveRDS(MOCK_DATA, paste0(OUTER_FOLDER, file_name))

#****************************************************************
# CANDADIAN DATA
#****************************************************************
data_type = 'Canadian'; seed_count = 1;

canadaX = c(1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 2, 3, 0, 2, 2,
            5, 7, 9, 7, 3, 4, 1, 4, 5, 7, 7, 7, 7, 3, 3, 5, 3, 5, 7, 4, 4, 2,
            3, 1, 1, 1, 0, 0, 2, 1, 3, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
            1, 1, 0, 1, 0, 0, 0, 1, 0, 2, 0, 3, 2, 2, 1, 2, 3, 4, 5, 5, 4, 6,
            6, 4, 8, 5, 6, 7, 5, 9, 1, 2, 3, 1, 1, 2, 0, 0, 0, 2, 0, 0, 0, 1)

#************************
#DATA I Extreme - No SS
canada_ss = rep(0, length(canadaX))
sim_data_canadaX1 = list(canadaX, canada_ss)

#************************
#DATA II Extreme SS; [1 0]
canada_bool = canadaX > 1
canada_ss = as.integer(canada_bool)
canada_ns = canadaX - canada_ss
sim_data_canadaX2 = list(canada_ns, canada_ss)

#PLOT
plot.ts(canadaX,
        main = 'SARs (2003) Canadian Outbreak',
        ylab = 'infection count')


#MODEL EVIDENCE DATA
#***********************
# 1. DATA 
#**********************
EPI_DATA = MOCK_DATA_6_DAYS
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
plot.ts(data_baseline, main = 'Baseline data')

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
EPI_DATA = data_baseline
EPI_DATA = MOCK_DATA_3_DAYS


#*************************
#SIMULATED DATA: LOAD
#*************************

# BASELINE DATA
OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
OUTER_FOLDER = paste0(OUTER_FOLDER, 'BASE_DATA/')
file_name = 'epi_data_base_1.rds'
data_baseline = readRDS(file = paste0(OUTER_FOLDER, file_name))

#SSEB DATA
file_name =  'epi_data_sseb_1.rds'
data_sseb = readRDS(file = paste0(DATA_FOLDER, file_name))
plot.ts(data_sseb)