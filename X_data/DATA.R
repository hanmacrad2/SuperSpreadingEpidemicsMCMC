#****************************************************************
# DATA
#****************************************************************
par(mfrow=c(2,1))
#library(SuperSpreadingEpidemicsMCMC)
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/X_data/DATA_SETS/"
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/X_data/DATA_SETS/DATA_SETS_50/"
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/X_data/"

#PARAMS
num_days = 50

#BASELINE
file_name = 'data_base_50_1_2.rds'
R0X = 1.2
data_baseline = SIMULATE_EPI_BASELINE(R0X, num_days = 50)
plot.ts(data_baseline, main = paste0('Baseline Data, R0 = ', R0X))
saveRDS(data_baseline, paste0(DATA_FOLDER, file_name))
data_baseline = readRDS(paste0(DATA_FOLDER, file_name))

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

#SSE
R0X = 1.2
kX = 0.1
file_name = 'data_sse_50_12_03.rds'
data_sse = SIMULATE_EPI_SSNB(num_days = num_days, R0 = R0X, k = kX)
plot.ts(data_sse, main = paste0('SSE Data, R0 = ', R0X, ', k = ', kX))

saveRDS(data_sse, paste0(DATA_FOLDER, file_name))

file_name = 'data_ssnb_1_6_2.rds'
data_ssnb2 = SIMULATE_EPI_SSNB()
plot.ts(EPI_DATA, main = 'SSE Data, R0 = 1.6')
saveRDS(data_ssnb, paste0(DATA_FOLDER, file_name))
data_ssnb = readRDS(paste0(DATA_FOLDER, file_name))

file_name = 'data_sse_1_3_1.rds'
data_sse = SIMULATE_EPI_SSNB(num_days = 30, R0 = 1.3, k = 1.0)
plot.ts(data_sse, main = 'SSE Data, R0 = 1.3, k = 1.0')
saveRDS(data_sse, paste0(DATA_FOLDER, file_name))

d2 = SIMULATE_EPI_SSNB(R0 = 1.2, k = 0.6)
plot.ts(d2)
                             
#SSI
R0X = 1.2
k = 0.1
file_name = 'data_ssi_50_1_2_0_3_v2.rds'
data_ssi = SIMULATE_EPI_SSIR(num_days = 50, R0X = R0X, k = k)
plot.ts(data_ssi$epidemic_data,
        main = paste0('SSI Data, R0: ', R0X, ', k: ', k),
        ylab = 'Daily count')
saveRDS(data_ssi, paste0(DATA_FOLDER, file_name))
data_ssi = readRDS(paste0(DATA_FOLDER, file_name))
EPI_DATA = data_ssi$epidemic_data

file_name = 'data_ssi.rds' 
data_ssir = SIMULATE_EPI_SSIR() #R0 - 1.6, k = 0.16
plot.ts(data_ssi$epidemic_data, main = 'SSI Data', ylab = 'Daily count')
saveRDS(data_ssir, paste0(DATA_FOLDER, file_name))
data_ssi = readRDS(paste0(DATA_FOLDER, file_name))


#SSE-B
file_name = 'data_sseb.rds'
data_sseb = SIMULATE_EPI_SSEB(num_days = 50)
plot.ts(data_sseb, main = 'SSE-B Data al: 0.8, be:0.2, ga:10' )#, R0 = 1.6')
EPI_DATA = data_sseb
saveRDS(data_sseb, paste0(DATA_FOLDER, file_name))
data_sseb = readRDS(paste0(DATA_FOLDER, file_name))

file_name = 'data_sseb_1_3.rds'
data_sseb = SIMULATE_EPI_SSEB(num_days = 50, betaX = 0.05)
plot.ts(data_sseb, main = 'SSE-B Data al: 0.8, be:0.05, ga:10' )#, R0 = 1.6')
EPI_DATA = data_sseb
saveRDS(data_sseb, paste0(DATA_FOLDER, file_name))
data_sseb = readRDS(paste0(DATA_FOLDER, file_name))

#SSIB
file_name = 'data_ssib.rds'
data_ssib = SIMULATE_EPI_SSIB(num_days = 50, aX = 0.8)
plot.ts(data_ssib,  main = 'SSIB DATA, R0 = 1.8', ylab = 'Daily infection count')
saveRDS(data_ssib, paste0(DATA_FOLDER, file_name))

file_name = 'data_ssib_50_1_3.rds'
data_ssib = SIMULATE_EPI_SSIB(num_days = 50, aX = 0.8, bX = 0.05)
plot.ts(data_ssib,  main = 'SSIB DATA, R0 = 1.3', ylab = 'Daily infection count')
saveRDS(data_ssib, paste0(DATA_FOLDER, file_name))
EPI_DATA = data_ssib

data_ssib = readRDS(paste0(DATA_FOLDER, file_name))
plot.ts(data_ssib)
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