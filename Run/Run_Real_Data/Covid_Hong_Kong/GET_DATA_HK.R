#*******************************
#1. DATA HONG KONG
#******************************
DATA_FOLDER = "~/Github/computing/REAL_DATA/1_HONG_KONG/DATA/"

#WAVE 1
file_name = 'df_hk_wave1.rds'
df_hk_wave1 = readRDS(file = paste0(DATA_FOLDER, file_name))
epi_hk1 = df_hk_wave1$cases
plot.ts(epi_hk1)
sum(epi_hk1)

#WAVE 2
file_name = 'df_hk_wave2.rds'
df_hk_wave2 = readRDS(file = paste0(DATA_FOLDER, file_name))
epi_hk2 = df_hk_wave2$cases
epidemic_data = epi_hk2
plot.ts(epi_hk2)
sum(epi_hk2)

#WAVE 3
file_name = 'df_hk_wave3.rds'
df_hk_wave3 = readRDS(file = paste0(DATA_FOLDER, file_name))
epi_hk3 = df_hk_wave3$cases
plot.ts(epi_hk3)
sum(epi_hk3)

#HK DATA TOTAL
file_name = 'df_hk.rds'
df_hk = readRDS(file = paste0(DATA_FOLDER, file_name))

#Plot
plot.ts(df_hk$cases)

#****************************************************************************
#* DATA PROCESS
#****************************************************************************

#****************************************************************************
#* WAVE 1
#****************************************************************************
end_date = "2020-04-03"
df_hk_wave1 <- df_hong_kong[df_hong_kong$date <= as.Date(end_date),]
plot.ts(df_hk_wave1$cases)

title = 'hong_kong Outbreak Canada, 2003. Wave 1.'
data_type = 'wave1_hong_kong_canada'
PLOT_EPI_DATA_DATE_PDF(df_hong_kong_wave1, RESULTS_FOLDER, title, data_type) 

#SAVE
file_name = 'df_hk_wave1.rds'
saveRDS(df_hk_wave1, file = paste0(DATA_FOLDER, file_name))
df_hk_wave1 = readRDS(file = paste0(DATA_FOLDER, file_name))

#EPI DATA
epi_hk1 = df_hk_wave1$cases
plot.ts(epi_hk1)

file_name = 'epi_hk1.rds'
saveRDS(epi_hk1, file = paste0(DATA_FOLDER, file_name))

#****************************************************************************
#* WAVE 2
#****************************************************************************

start_date <- as.Date("2020-06-15")
end_date <- as.Date("2020-08-01") 

#start_date <- as.Date("2020-05-17")

# Using subset function
df_hk_wave2 <- subset(df_hk, date >= start_date & date <= end_date)
plot.ts(df_hk_wave2$cases)
epi_hk2 = df_hk_wave2$cases
plot.ts(epi_hk2)

#SAVE
file_name = 'df_hk_wave_2.rds'
saveRDS(df_hk_wave2, file = paste0(DATA_FOLDER, file_name))

file_name = 'epi_hk2.rds'
saveRDS(epi_hk2, file = paste0(DATA_FOLDER, file_name))

#PLOT
title = 'hong_kong Outbreak Canada, 2003. Wave 2.'
data_type = 'wave2_hong_kong_canada'
PLOT_EPI_DATA_DATE_PDF(df_hong_kong_wave1, RESULTS_FOLDER, title, data_type) 

#******************
#* WAVE 3
start_date <- as.Date("2020-10-28")
end_date <- as.Date("2020-12-19")

df_hk_wave3 <- subset(df_hk, date >= start_date & date <= end_date)
plot.ts(df_hk_wave3$cases)

#SAVE
file_name = 'df_hk_wave3.rds'
saveRDS(df_hk_wave3, file = paste0(DATA_FOLDER, file_name))
df_hk_wave3 = readRDS(file = paste0(DATA_FOLDER, file_name))

#EPI DATA
epi_hk3 = df_hk_wave3$cases
plot.ts(epi_hk3)

file_name = 'epi_hk3.rds'
saveRDS(epi_hk3, file = paste0(DATA_FOLDER, file_name))


#******************
#* WAVE 4 (WAVE 3 EDITED)
start_date <- as.Date("2020-11-07")
end_date <- as.Date("2020-12-19")

df_hk_wave4 <- subset(df_hk, date >= start_date & date <= end_date)
plot.ts(df_hk_wave4$cases)

#SAVE
file_name = 'df_hk_wave4.rds'
saveRDS(df_hk_wave4, file = paste0(DATA_FOLDER, file_name))
df_hk_wave4 = readRDS(file = paste0(DATA_FOLDER, file_name))

#EPI DATA
epi_hk4 = df_hk_wave4$cases
plot.ts(epi_hk4)

file_name = 'epi_hk4.rds'
saveRDS(epi_hk4, file = paste0(DATA_FOLDER, file_name))

file_name = 'epi_hk3.rds'
saveRDS(epi_hk3, file = paste0(DATA_FOLDER, file_name))


