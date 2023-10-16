#PLOT INFERENCE RESULTS

#******************
#* 1. BASELINE 
PLOT_CI_PARAMS(df_base, COMP_FOLDER, PDF = FALSE, 
               FLAG_PARAM = list(r0 = TRUE, k = FALSE, 
                                 alpha = FALSE, gamma = FALSE),
               FLAG_MODEL = list(BASELINE = TRUE, SSE = FALSE, SSI = TRUE,
                                 SSEB = FALSE, SSIB = FALSE),
               FLAG_FILTER = list(end_day = FALSE, tot_infs = TRUE))

#******************
#* 1. SSE
PLOT_CI_PARAMS(df_sse, COMP_FOLDER, PDF = FALSE, 
               FLAG_PARAM = list(r0 = FALSE, k = TRUE, 
                                 alpha = FALSE, gamma = FALSE),
               FLAG_MODEL = list(BASELINE = FALSE, SSE = TRUE, SSI = FALSE,
                                 SSEB = FALSE, SSIB = FALSE),
               FLAG_FILTER = list(end_day = FALSE, tot_infs = TRUE))

#****************
#FIXED
PLOT_CI_PARAMS_FIXED(df_sse1, COMP_FOLDER, fig_num = '01', PDF = TRUE, 
                                 FIXED_PARAM = list(r0 = TRUE, k = FALSE,
                                                    alpha = FALSE, gamma = FALSE),
                                 VAR_PARAM = list(r0 = FALSE, k = TRUE,
                                                  alpha = FALSE, gamma = FALSE),
                                 FLAG_MODEL = list(BASELINE = FALSE, SSE = TRUE, SSI = FALSE,
                                                   SSEB = FALSE, SSIB = FALSE))


#******************
#* 1. SSI 
PLOT_CI_PARAMS(df_ssi1, COMP_FOLDER, PDF = FALSE, 
               FLAG_PARAM = list(r0 = FALSE, k = TRUE, 
                                 alpha = FALSE, gamma = FALSE),
               FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = TRUE,
                                 SSEB = FALSE, SSIB = FALSE),
               FLAG_FILTER = list(end_day = FALSE, tot_infs = TRUE))
