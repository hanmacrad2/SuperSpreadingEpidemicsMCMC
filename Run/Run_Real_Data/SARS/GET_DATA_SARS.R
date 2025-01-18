#REAL DATA
library(outbreaks)
library(dplyr)

#*******************************
#1. SARS CANADA 
#******************************
#DATA
DATA_FOLDER = "~/Github/computing/REAL_DATA/2_SARS/DATA/"

data("sars_canada_2003")
df_sars = sars_canada_2003

df_sars <- df_sars %>%
  mutate(cases = cases_travel + cases_household + cases_healthcare + cases_other)

epi_data_sars = df_sars$total_cases
plot.ts(epi_data_sars)

#Save
file_name = 'df_sars.rds'
saveRDS(df_sars, file = paste0(DATA_FOLDER, file_name))

df_sars = readRDS(file = paste0(DATA_FOLDER, file_name))
  
#SAVE
file_name = 'epi_data_sars_canada_03.rds'
saveRDS(epi_data_sars, file = paste0(DATA_FOLDER, file_name))

#************************************
#* WAVE 1
#************************************
end_date = "2003-04-04"
df_sars_wave1 <- df_sars[df_sars$date <= as.Date(end_date),]

title = 'SARS Outbreak Canada, 2003. Wave 1.'
data_type = 'wave1_sars_canada'
PLOT_EPI_DATA_DATE_PDF(df_sars_wave1, RESULTS_FOLDER, title, data_type) 

#SAVE
file_name = 'df_sars_wave1.rds'
saveRDS(df_sars_wave1, file = paste0(DATA_FOLDER, file_name))

df_sars_wave1 = readRDS(file = paste0(DATA_FOLDER, file_name)) 

#EPI DATA
epi_wave1 = df_sars_wave1$cases
plot.ts(epi_wave1)
epidemic_data = epi_wave1

file_name = 'epi_wave1.rds'
saveRDS(epi_wave1, file = paste0(DATA_FOLDER, file_name))

epi_wave1_sars = readRDS(file = paste0(DATA_FOLDER, file_name)) 
plot.ts(epi_wave1_sars)

#************************************
#* WAVE 2
#************************************
start_date = "2003-05-09"
end_date = "2003-05-29"
df_sars_wave2 <- subset(df_sars, date >= start_date & date <= end_date)

title = 'Wave 2. SARS Outbreak Canada, 2003'
data_type = 'wave2_sars_canada'
PLOT_EPI_DATA_DATE_PDF(df_sars_wave2, RESULTS_FOLDER, title, data_type) 

#SAVE
DATA_FOLDER = paste0(RESULTS_FOLDER, '/DATA/')
create_folder(DATA_FOLDER)
file_name = 'df_sars_wave2.rds'
saveRDS(df_sars_wave2, file = paste0(DATA_FOLDER, file_name))

file_name = 'df_sars_wave2.rds'
df_sars_wave2 = readRDS(file = paste0(DATA_FOLDER, file_name))

#EPI DATA
epi_wave2 = df_sars_wave2$cases
plot.ts(epi_wave2)

epidemic_data = epi_wave2

file_name = 'epi_wave2.rds'
saveRDS(epi_wave2, file = paste0(DATA_FOLDER, file_name))
epi_wave2 = readRDS(file = paste0(DATA_FOLDER, file_name))
plot.ts(epi_wave2)
