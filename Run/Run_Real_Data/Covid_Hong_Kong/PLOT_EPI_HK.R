#***********************************************************************************
#* 1. PLOT EPI DATA HONG KONG 
#********************************************************************************
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/1_HONG_KONG/'

#* LOAD_EPI DATA
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/HONG_KONG/DATA/'
file_name = 'df_hk.rds'
df_hk = readRDS(file = paste0(DATA_FOLDER, file_name))

#PLOT 
RESULTS_FOLDER = "~/Github/computing/REAL_DATA/HONG_KONG/plots/"
title = 'SARS-CoV-2 Outbreak Hong Kong, 2020-2022. Daily Cases'
data_type = 'Hong_kong_'
PLOT_EPI_DATA_DATE_PDF(df_hk, RESULTS_FOLDER, title, data_type) 

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

#*****************************************************************
#* PLOT WAVES 
#****************************************************************

#PLOT
RESULTS_FOLDER = "~/Github/computing/REAL_DATA/1_HONG_KONG/RESULTS/"
data_type = 'waves_hk'
GET_EPI_DATA_PDF(RESULTS_FOLDER, data_type)
par(mfrow=c(3,1))
par(oma = c(1, 1, 1, 1))
par(mar = c(4.5,5,4,4))
lwd_data = 1.7

ylim = c(0, max(df_hk_wave1$cases))
title = 'Wave 1. Daily Cases: January 23rd - April 4th 2020' #April 26th
PLOT_EPI_DATA_DATE(df_hk_wave1, title, cex = 1.5, ylim=ylim, lwd_data = lwd_data)

ylim = c(0, max(df_hk_wave2$cases))
title = 'Wave 2. Daily Cases: June 15th - August 1st 2020'
PLOT_EPI_DATA_DATE(df_hk_wave2, title, cex = 1.5, ylim=ylim, lwd_data = lwd_data)

ylim = c(0, max(df_hk_wave3$cases))
title = 'Wave 3. Daily Cases: October 28th - December 19th 2020' 
PLOT_EPI_DATA_DATE(df_hk_wave3, title, cex = 1.5, ylim=ylim,lwd_data = lwd_data)

dev.off()



#************************************
#* EPI DATA
#***********************************
#PLOT
par(mfrow=c(3,1))
#cex = 1.0; main_font = 2.5; axis_font = 1.6

title = 'Wave 1. January 23rd - April 26th 2020'
PLOT_EPI_DATA(epidemic_data_wave1, title, cex = 1.5)

title = 'Wave 2. May 17th - August 1st 2020'
PLOT_EPI_DATA(epidemic_data_wave2, title, cex = 1.5)

title = 'Wave 3. September 23rd - December 19th 2020'
PLOT_EPI_DATA(epidemic_data_wave3, title, cex = 1.5)

dev.off()