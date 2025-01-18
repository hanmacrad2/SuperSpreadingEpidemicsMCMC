#***********************************************************
#* DATA NZ
#**********************************************************

DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'

#DATA 1 NZ SOUTH, 2020
file_name = 'df_nz_south_20.rds'
df_nz_south = readRDS(file= paste0(DATA_FOLDER, file_name))
plot(df_nz_south$dates, df_nz_south$cases, type = 'l')

data_type = '  SARS-CoV-2 - New Zealand South Outbreak, 2020'
DATA_SET = 'SARS-CoV-2, NZ South, 2020'

#DATA 2 NZ - WAITEMATA, AUCKLAND, 2021
file_name = 'df_nz_waita_21.rds'
df_nz_waita = readRDS(file= paste0(DATA_FOLDER, file_name))
plot(df_nz_waita$dates, df_nz_waita$cases, type = 'l')

epidemic_data = df_nz_south$cases
epidemic_data = df_nz_waita$cases
data_type = 'SARS CoV-2 Outbreak - Auckland district, New Zealand 2021'
DATA_SET = 'SARS-CoV-2 - Auckland NZ, 2021'

#**************************************************************
#DATA NZ TOTAL
#**************************************************************
library(readxl)
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'

#DATA TOTAL
file_name = 'nz_covid_cases_counts_location.xlsx'
df_nz_total = read_excel(paste0(DATA_FOLDER, file_name))

#FILTER
names(df_nz_total)[names(df_nz_total) == "Report date"] <- "date"
names(df_nz_total)[names(df_nz_total) == "Cases"] <- "cases"
df_nz_total <- df_nz_total[df_nz_total$date != "Total", ]
df_nz_total$date <- as.Date(df_nz_total$date, format = "%Y-%m-%d")

#SAVE
file_name = 'df_nz_total.rds'
saveRDS(df_nz_total, file = paste0(DATA_FOLDER, file_name))
df_nz_total = readRDS(file = paste0(DATA_FOLDER, file_name))

#****************************************************************************
#* DATA: 2020 & 2021
#****************************************************************************
end_date = "2021-12-30"
df_nz_20 <- df_nz_total[df_nz_total$date <= as.Date(end_date),]
plot(df_nz_20$date, df_nz_20$cases, type = 'l')

file_name = 'df_nz_20.rds'
saveRDS(df_nz_20, file = paste0(DATA_FOLDER, file_name))

#PLOT
lwd_data = 1.7
title = 'SARS CoV-2, New Zealand Daily Cases. 2020-2021'
data_type = 'nz_total_data_20_21'

PLOT_EPI_DATA_DATE_PDF(df_nz_20, RESULTS_FOLDER, title, data_type,
                       lwd_data = lwd_data, X_AXIS_DATES = TRUE)

#**********************************************************************************
# DATA NZ WEDDING 2020
#********************************************************************************

#DATA RUN 2
dates <- as.Date(c( "2020-03-26", 
                    "2020-03-27", "2020-03-28", "2020-03-29", "2020-03-30", "2020-03-31", 
                    "2020-04-01", "2020-04-02", "2020-04-03", "2020-04-04", "2020-04-05", 
                    "2020-04-06"))
cases <- c(2, 16, 8, 12, 7, 7, 12, 15, 11, 8, 16, 8)

# Create the dataframe
df_nz <- data.frame(dates = dates, cases = cases)
df_nz_south = df_nz

#SAVE
file_name = 'df_nz_south_20.rds'
saveRDS(df_nz_south, file= paste0(DATA_FOLDER, file_name))
df_nz = readRDS(file= paste0(DATA_FOLDER, file_name))

#EPIDEMIC DATA
epidemic_data = df_nz$cases
plot.ts(epidemic_data)
file_name = 'epi_data_nz2.rds'
saveRDS(epidemic_data, file= paste0(DATA_FOLDER, file_name))
epidemic_data = readRDS(file= paste0(DATA_FOLDER, file_name))


#*****************************************************************
#* 
#*  DATA WAITEMATA (DATA SET 2)
#*  
#* ******************************************************************
#DATA RUN 2
#DATES: 
#2021-08-17: NZ moves to Alert level 4
#2021-08-23: Only Auckland remains at alert level 4. Alert level 4 being the highest level of restrictions

cases =  c(5, 3, 9, 4, 5, 11, 7, 22, 5, 14, 13, 14, 17, 14, 20, 4)

# Create vectors for the dates and cases
dates <- as.Date(c("2021-08-17", "2021-08-18", "2021-08-19", "2021-08-20", "2021-08-21", 
                   "2021-08-22", "2021-08-23", "2021-08-24", "2021-08-25", "2021-08-26", 
                   "2021-08-27", "2021-08-28", "2021-08-29", "2021-08-30", "2021-08-31", 
                   "2021-09-01"))

# Create a dataframe
df_nz_waita2 <- data.frame(dates, cases)

#SAVE
file_name = 'df_nz_waita2.rds'
saveRDS(df_nz_waita2, file= paste0(DATA_FOLDER, file_name))
df_nz_waita2 = readRDS(file= paste0(DATA_FOLDER, file_name))

epidemic_data = df_nz_waita2$cases
plot.ts(epidemic_data)
file_name = 'epi_data_waita2.rds'
saveRDS(epidemic_data, file= paste0(DATA_FOLDER, file_name))

#******************************
#DATA -> PLOT
# Create vectors for the dates and cases
dates_plot <- as.Date(c("2021-08-14", "2021-08-15", "2021-08-16", "2021-08-17", "2021-08-18", "2021-08-19", "2021-08-20", "2021-08-21", 
                        "2021-08-22", "2021-08-23", "2021-08-24", "2021-08-25", "2021-08-26", 
                        "2021-08-27", "2021-08-28", "2021-08-29", "2021-08-30", "2021-08-31", 
                        "2021-09-01", "2021-09-02", "2021-09-03", "2021-09-04", "2021-09-05", 
                        "2021-09-06", "2021-09-07", "2021-09-08"))

cases_plot <- c(0, 0, 0, 5, 3, 9, 4, 5, 11, 7, 22, 5, 14, 13, 14, 17, 14, 20, 4, 2, 10, 5, 2, 4, 3, 8)

# Create a dataframe
df_nz_waita_plot <- data.frame(dates_plot, cases_plot)

#PLOT
title = 'SARS CoV-2 Outbreak - Auckland district, New Zealand 2021'
data_type = 'nz_sars_cov2_waiteama_2021'
PLOT_EPI_DATA_DATE_PDF(df_nz_waita_plot, RESULTS_FOLDER, title, data_type, SSE = TRUE) 


