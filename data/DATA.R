#****************************************************************
# DATA
#****************************************************************
#library(SuperSpreadingEpidemicsMCMC)
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"

#************************
#SIMULATE DATA
#************************

#BASELINE
file_name = 'data_baseline.rds'
data_baseline = readRDS(paste0(DATA_FOLDER, file_name))
plot.ts(data_baseline)

data_baseline4 = SIMULATE_EPI_BASELINE(2.1, num_days = 30)
plot.ts(data_baseline4)
saveRDS(data_baseline2, paste0(DATA_FOLDER, file_name))

#SSNB
file_name = 'data_ssnb.rds'
data_ssnb = SIMULATE_EPI_SSNB(num_days = 50)
plot.ts(data_ssnb)
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
MOCK_DATA_3_DAYS = c(1,2,1)
MOCK_DATA_3_DAYS = c(1,1,1)
MOCK_DATA_4_DAYS = c(1,1,0,1)
MOCK_DATA_6_DAYS = c(1,1,0,1,0,1)
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