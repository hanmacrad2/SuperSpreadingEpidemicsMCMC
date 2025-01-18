#SARS BOTH WAVES
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/2_SARS/RESULTS/'
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/2_SARS/DATA/'

#DATA
file_name = 'df_sars_wave1.rds'
df_sars_wave1 = readRDS(paste0(DATA_FOLDER, file_name))

file_name = 'df_sars_wave2.rds'
df_sars_wave2 = readRDS(paste0(DATA_FOLDER, file_name))
#epidemic_data2 = df_sars_wave2$cases
#plot.ts(epidemic_data2)

#**********************************************************
#* #PLOT TWO WAVES TOGETHER: EPI DATA + MODEL SELECTION RESULTS 


#SETTING
cex = 1.8
data_type = 'waves_both_sars_canada'
result_type = 'epi_data_model_comp'

GET_PDF_SETTING(RESULTS_FOLDER, data_type, result_type, plot_width = 17.0,  
                plot_height = 11.0)

par(mfrow=c(2,2))
par(oma = c(1, 1, 1, 1))
par(mar = c(4.5,5,4,4))

#EPI DATA
#WAVE 1
title = 'Wave 1. SARS Outbreak Canada, 2003'
ylim = c(0, max(df_sars_wave1$cases))
PLOT_EPI_DATA_DATE(df_sars_wave1, title, cex = cex, ylim = ylim)

#WAVE 2
title = 'Wave 2. SARS Outbreak Canada, 2003'
ylim = c(0, max(df_sars_wave2$cases))
PLOT_EPI_DATA_DATE(df_sars_wave2, title, cex = cex, ylim = ylim)

#***********************
#PLOT POST PROBS
post_probs1 = c(0,0, 0.23 ,0, 0.77)
post_probs2 = c(0,0, 0.20 ,0, 0.80)

title = 'Wave 1 SARS, 2003.'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs1[1], SSE = post_probs1[2],
                         SSI = post_probs1[3], SSEB = post_probs1[4], SSIB = post_probs1[5]),
                    title = title, cex = cex)

#PLOT POST PROBS
title = 'Wave 2 SARS, 2003.'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs2[1], SSE = post_probs2[2],
                         SSI = post_probs2[3], SSEB = post_probs2[4], SSIB = post_probs2[5]),
                    title = title, cex = cex)

dev.off()

