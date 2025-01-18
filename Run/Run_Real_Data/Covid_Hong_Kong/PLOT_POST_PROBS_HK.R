#SARS BOTH WAVES
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/1_HONG_KONG/RESULTS/'
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/1_HONG_KONG/DATA/'
#DATA
#df_hk_wave1, df_hk_wave2, df_hk_wave4

#**********************************************************
#* #PLOT THREE WAVES TOGETHER: EPI DATA + MODEL SELECTION RESULTS 
lwd_data = 2.25

#SETTING
cex = 1.8
data_type = 'hk_waves'
result_type = 'epi_data_model_comp'

GET_PDF_SETTING(RESULTS_FOLDER, data_type, result_type, plot_width = 17.0,  
                plot_height = 8.5)

par(mfrow=c(2,3))
par(oma = c(1, 1, 1, 1))
par(mar = c(4.5,5,4,4))

#EPI DATA
#WAVE 1
ylim = c(0, max(df_hk_wave1$cases))
#title = 'Wave 1. Daily Cases Hon'
title = 'SARS-CoV-2 Outbreak, Wave 1 Hong Kong 2020'
PLOT_EPI_DATA_DATE(df_hk_wave1, title, cex = cex, ylim = ylim, lwd_data = lwd_data)

#WAVE 2
ylim = c(0, max(df_hk_wave2$cases))
title = 'SARS-CoV-2 Outbreak, Wave 2 Hong Kong 2020'
PLOT_EPI_DATA_DATE(df_hk_wave2, title, cex = cex, ylim = ylim, lwd_data = lwd_data)

#WAVE 3
ylim = c(0, max(df_hk_wave4$cases))
title = 'SARS-CoV-2 Outbreak, Wave 3 Hong Kong 2020'
PLOT_EPI_DATA_DATE(df_hk_wave4, title, cex = cex, ylim = ylim, lwd_data = lwd_data)

#*************************
#PLOT POST PROBS
DATA_SET = 'SARS-CoV-2, Wave 1 Hong Kong, 2020'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs[1], SSE = post_probs[2],
                         SSI = post_probs[3], SSEB = post_probs[4], SSIB = post_probs[5]),
                    title = "", cex = cex) # title = DATA_SET

#PLOT POST PROBS
DATA_SET = 'SARS-CoV-2, Wave 2 Hong Kong, 2020'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs2[1], SSE = post_probs2[2],
                         SSI = post_probs2[3], SSEB = post_probs2[4], SSIB = post_probs2[5]),
                    title = "", cex = cex)

#PLOT POST PROBS
DATA_SET = 'SARS-CoV-2, Wave 3 Hong Kong, 2020'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs3[1], SSE = post_probs3[2],
                         SSI = post_probs3[3], SSEB = post_probs3[4], SSIB = post_probs3[5]),
                    title = "", cex = cex)

dev.off()


#POST PROBS
post_probs = c(0, 0.9995, 0, 0.0005, 0)
post_probs2 = c(0, 0.97, 0, 0.03, 0)
post_probs3 = c(0, 0.996, 0.001, 0.003, 0)
#post_probs3 = c(0, 0.996, 0.003, 0.001, 0)
