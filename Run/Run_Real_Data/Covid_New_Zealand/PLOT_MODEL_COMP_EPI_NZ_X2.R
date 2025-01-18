#SARS BOTH WAVES
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/RESULTS/'
DATA_FOLDER = '~/GitHub/computing/REAL_DATA/3_NZ/DATA/'

#DATA
#file_name = 'df_nz.rds'
file_name = 'df_nz_wedding_plot_model_comp.rds'
df_nz = readRDS(paste0(DATA_FOLDER, file_name))

#file_name = 'df_nz_waita2.rds'
file_name = 'df_nz_waita_plot_model_comp.rds'
df_nz_waita = readRDS(paste0(DATA_FOLDER, file_name))

#**********************************************************
#* #PLOT TWO WAVES TOGETHER: EPI DATA + MODEL SELECTION RESULTS 

post_probs1 = c(0, 0.68, 0, 0.32, 0)
post_probs2 = c(0, 0.75, 0, 0.25, 0)

#SETTING
cex = 1.8
data_type = 'nz_both_datasets'
result_type = 'epi_data_model_comp'

GET_PDF_SETTING(RESULTS_FOLDER, data_type, result_type, plot_width = 17.0,  
                plot_height = 11.0)

par(mfrow=c(2,2))
par(oma = c(1, 1, 1, 1))
par(mar = c(4.5,5,4,4))

#EPI DATA
#WAVE 1
ylim = c(0, max(df_nz$cases))
title = 'SARS-CoV-2, Outbreak New Zealand South, 2020'
PLOT_EPI_DATA_DATE(df_nz, title, cex = cex, ylim = ylim)

#WAVE 2
ylim = c(0, max(df_nz_waita$cases))
title = 'SARS CoV-2 Outbreak, Auckland New Zealand, 2021'
PLOT_EPI_DATA_DATE(df_nz_waita, title, cex = cex, ylim = ylim)

#PLOT POST PROBS
DATA_SET = 'SARS-CoV-2, NZ South, 2020'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs1[1], SSE = post_probs1[2],
                         SSI = post_probs1[3], SSEB = post_probs1[4], SSIB = post_probs1[5]),
                    title = "", cex = cex) # title = DATA_SET

#PLOT POST PROBS
DATA_SET = 'SARS-CoV-2 - Auckland NZ, 2021'
BAR_PLOT_POST_PROBS(list(Baseline = post_probs2[1], SSE = post_probs2[2],
                         SSI = post_probs2[3], SSEB = post_probs2[4], SSIB = post_probs2[5]),
                    title = "", cex = cex)

dev.off()

