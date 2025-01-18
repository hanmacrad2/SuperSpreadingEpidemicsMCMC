#***********************************************************************************
#* 1. PLOT EPI DATA HONG KONG 
#********************************************************************************
#DATA NZ TOTAL
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'
file_name = 'df_nz_20.rds'
df_nz_20 = readRDS(file = paste0(DATA_FOLDER, file_name))

#PLOT
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/'
title = 'SARS CoV-2, New Zealand Daily Cases. 2020-2021'
data_type = 'nz_total_data_20_21'

PLOT_EPI_DATA_DATE_PDF(df_nz_20, RESULTS_FOLDER, title, data_type,
                       lwd_data = lwd_data, plot_width = 15.5,  
                       plot_height = 7.0, cex = 1.6,
                       X_AXIS_DATES = TRUE)

#************************************************************
#* WAVES 
#***********************************************************
file_name = 'df_hk_wave1.rds'
df_hk_wave1 = readRDS(file = paste0(DATA_FOLDER, file_name))

file_name = 'df_hk_wave2.rds'
df_hk_wave2 = readRDS(file = paste0(DATA_FOLDER, file_name))

file_name = 'df_hk_wave3.rds'
file_name = 'df_hk_wave3_short.rds'
df_hk_wave3 = readRDS(file = paste0(DATA_FOLDER, file_name))


#*****************************
#* DATA
#* ****************************


#******************************
#DATA -> PLOT

#1. NZ WEDDING
dates_plot <- as.Date(c("2020-03-18", "2020-03-19", "2020-03-20", "2020-03-21", "2020-03-22", "2020-03-23", "2020-03-24", "2020-03-25", "2020-03-26", 
                        "2020-03-27", "2020-03-28", "2020-03-29", "2020-03-30", "2020-03-31", 
                        "2020-04-01", "2020-04-02", "2020-04-03", "2020-04-04", "2020-04-05", 
                        "2020-04-06"))
cases_plot <- c(0,0, 0, 0, 1, 3, 2, 5, 2, 16, 8, 12, 7, 7, 12, 15, 11, 8, 16, 8)

# Create the dataframe
df_nz_plot <- data.frame(date = dates_plot, cases = cases_plot)

#2. NZ WAITA
dates_plot <- as.Date(c("2021-08-14", "2021-08-15", "2021-08-16", "2021-08-17", "2021-08-18", "2021-08-19", "2021-08-20", "2021-08-21", 
                        "2021-08-22", "2021-08-23", "2021-08-24", "2021-08-25", "2021-08-26", 
                        "2021-08-27", "2021-08-28", "2021-08-29", "2021-08-30", "2021-08-31", 
                        "2021-09-01"))

cases_plot =  c(0, 0, 0, 5, 3, 9, 4, 5, 11, 7, 22, 5, 14, 13, 14, 17, 14, 20, 4)

# Create a dataframe
df_nz_waita_plot <- data.frame(dates_plot, cases_plot)

#*****************************************************************
#* PLOT WAVES 
#****************************************************************

#PLOT
RESULTS_FOLDER = "~/Github/computing/REAL_DATA/3_NZ/RESULTS/"
data_type = 'nz'
GET_EPI_DATA_PDF(RESULTS_FOLDER, data_type, plot_width = 11.5, plot_height = 8.2) 
par(mfrow=c(2,1))
par(oma = c(1, 1, 1, 1))
par(mar = c(4.5,5,4,4))
lwd_data = 1.7

ylim = c(0, max(df_nz_plot$cases))
title = 'SARS CoV-2 Outbreak & SSE (Wedding) - New Zealand South 2020'
PLOT_EPI_DATA_DATE(df_nz_plot, title, cex = 1.2, ylim=ylim, lwd_data = lwd_data, SSE = TRUE)

ylim = c(0, max(df_nz_waita_plot$cases))
title =  'SARS CoV-2 Outbreak - Auckland district, New Zealand 2021' 
PLOT_EPI_DATA_DATE(df_nz_waita2, title, cex = 1.2, ylim=ylim, lwd_data = lwd_data, SSE = TRUE)


dev.off()
