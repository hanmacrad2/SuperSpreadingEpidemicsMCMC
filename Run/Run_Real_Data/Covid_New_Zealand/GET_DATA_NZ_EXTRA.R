#NZ DATA

#************************************************************
#* DATA 1 - NZ WEDDING
#***********************************************************
#DATA
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'

#SAVE
file_name = 'df_nz.rds'
df_nz = readRDS(file= paste0(DATA_FOLDER, file_name))

#EPIDEMIC DATA
epidemic_data = df_nz$cases
plot.ts(epidemic_data)
file_name = 'epi_data_nz.rds'
saveRDS(epidemic_data, file= paste0(DATA_FOLDER, file_name))

#********************************************************************************
#DATA - TO PLOT
#********************************************************************************
dates_plot <- as.Date(c("2020-03-18", "2020-03-19", "2020-03-20", "2020-03-21", "2020-03-22", "2020-03-23", "2020-03-24", "2020-03-25", "2020-03-26", 
                   "2020-03-27", "2020-03-28", "2020-03-29", "2020-03-30", "2020-03-31", 
                   "2020-04-01", "2020-04-02", "2020-04-03", "2020-04-04", "2020-04-05", 
                   "2020-04-06"))
cases_plot <- c(0,0, 0, 0, 1, 3, 2, 5, 2, 16, 8, 12, 7, 7, 12, 15, 11, 8, 16, 8)

# Create the dataframe
df_nz_plot <- data.frame(date = dates_plot, cases = cases_plot)

#PLOT
title = 'SARS CoV-2 Outbreak & SSE (Wedding) - New Zealand South 2020'
data_type = 'nz_2020_sars_cov_2'
PLOT_EPI_DATA_DATE_PDF(df_nz_plot, RESULTS_FOLDER, title, data_type, SSE = TRUE) 

#********************************************************************************
#DATA - TO PLOT - MODEL COMPARISON !! 
#********************************************************************************
dates_plot <- as.Date(c("2020-03-21", "2020-03-22", "2020-03-23", "2020-03-24", "2020-03-25", "2020-03-26", 
                        "2020-03-27", "2020-03-28", "2020-03-29", "2020-03-30", "2020-03-31", 
                        "2020-04-01", "2020-04-02", "2020-04-03", "2020-04-04", "2020-04-05", 
                        "2020-04-06"))
cases_plot <- c(0, 1, 3, 2, 5, 2, 16, 8, 12, 7, 7, 12, 15, 11, 8, 16, 8)

# Create the dataframe
df_nz_plot_model_comp <- data.frame(date = dates_plot, cases = cases_plot)

#SAVE 
file_name = 'df_nz_wedding_plot_model_comp.rds'
saveRDS(df_nz_plot_model_comp, file = paste0(DATA_FOLDER, file_name))

#PLOT
title = 'SARS CoV-2 Outbreak & SSE (Wedding) - New Zealand South 2020'
data_type = 'nz_2020_sars_cov_2'
PLOT_EPI_DATA_DATE_PDF(df_nz_plot, RESULTS_FOLDER, title, data_type, SSE = TRUE) 

#**********************************************************************************
# DATA NZ USE 
#********************************************************************************

#DATA RUN 2
dates <- as.Date(c( "2020-03-26", 
                    "2020-03-27", "2020-03-28", "2020-03-29", "2020-03-30", "2020-03-31", 
                    "2020-04-01", "2020-04-02", "2020-04-03", "2020-04-04", "2020-04-05", 
                    "2020-04-06"))
cases <- c(2, 16, 8, 12, 7, 7, 12, 15, 11, 8, 16, 8)

# Create the dataframe
df_nz2 <- data.frame(date = dates, cases = cases)

#SAVE
file_name = 'df_nz2.rds'
saveRDS(df_nz2, file= paste0(DATA_FOLDER, file_name))

#EPIDEMIC DATA
epidemic_data = df_nz2$cases
plot.ts(epidemic_data)
file_name = 'epi_data_nz2.rds'
saveRDS(epidemic_data, file= paste0(DATA_FOLDER, file_name))
epidemic_data = readRDS(file= paste0(DATA_FOLDER, file_name))

#*****************************************************
#DATA RUN 1
#****************************************************
dates <- as.Date(c( "2020-03-23", "2020-03-24", "2020-03-25", "2020-03-26", 
                   "2020-03-27", "2020-03-28", "2020-03-29", "2020-03-30", "2020-03-31", 
                   "2020-04-01", "2020-04-02", "2020-04-03", "2020-04-04", "2020-04-05", 
                   "2020-04-06"))
cases <- c(3, 2, 5, 2, 16, 8, 12, 7, 7, 12, 15, 11, 8, 16, 8)
df_nz <- data.frame(date = dates, cases = cases)

#SAVE
file_name = 'df_nz.rds'
saveRDS(df_nz, file= paste0(DATA_FOLDER, file_name))

#EPIDEMIC DATA
epidemic_data = df_nz$cases
plot.ts(epidemic_data)
file_name = 'epi_data_nz.rds'
saveRDS(epidemic_data, file= paste0(DATA_FOLDER, file_name))

#*****************************************************************
#* 
#*  DATA WAITEMATA (DATA SET 2) USED - SEE GETA DATA NZ
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

#*************************************************************************
#*
#* DATA WAITEMATA RUN 1
#* 
#*********************************************************************
cases <- c(5, 3, 9, 4, 5, 11, 7, 22, 5, 14, 13, 14, 17, 14, 20, 4, 2, 10, 5, 2, 4, 3, 8)

dates_plot <- as.Date(c("2021-08-17", "2021-08-18", "2021-08-19", "2021-08-20", "2021-08-21", 
                        "2021-08-22", "2021-08-23", "2021-08-24", "2021-08-25", "2021-08-26", 
                        "2021-08-27", "2021-08-28", "2021-08-29", "2021-08-30", "2021-08-31", 
                        "2021-09-01", "2021-09-02", "2021-09-03", "2021-09-04", "2021-09-05", 
                        "2021-09-06", "2021-09-07", "2021-09-08"))

# Create a dataframe
df_nz_waita <- data.frame(dates, cases)

#SAVE
file_name = 'df_nz_waita.rds'
saveRDS(df_nz_waita, file= paste0(DATA_FOLDER, file_name))
df_nz_waita = readRDS(file= paste0(DATA_FOLDER, file_name))

epidemic_data = df_nz_waita$cases
plot.ts(epidemic_data)
file_name = 'epi_data_waita.rds'
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


#********************************************************************************
#DATA - TO PLOT - MODEL COMPARISON !! 
#********************************************************************************

cases_plot =  c(0, 5, 3, 9, 4, 5, 11, 7, 22, 5, 14, 13, 14, 17, 14, 20, 4)

# Create vectors for the dates and cases
dates_plot <- as.Date(c("2021-08-16", "2021-08-17", "2021-08-18", "2021-08-19", "2021-08-20", "2021-08-21", 
                   "2021-08-22", "2021-08-23", "2021-08-24", "2021-08-25", "2021-08-26", 
                   "2021-08-27", "2021-08-28", "2021-08-29", "2021-08-30", "2021-08-31", 
                   "2021-09-01"))

# Create the dataframe
df_nz_waita_plot_model_comp <- data.frame(date = dates_plot, cases = cases_plot)

#SAVE 
file_name = 'df_nz_waita_plot_model_comp.rds'
saveRDS(df_nz_waita_plot_model_comp, file = paste0(DATA_FOLDER, file_name))
