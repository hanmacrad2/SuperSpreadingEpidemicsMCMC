#PLOT INFERENCE RESULTS

#******************
#* 1. BASELINE
COMP_FOLDER = "~/GitHub/computing/mcmc/baseline/"
PLOT_CI_PARAMS(df_base, COMP_FOLDER, fig_num = '01', PDF = TRUE,  
               FLAG_PARAM = list(r0 = TRUE, k = FALSE, 
                                 alpha = FALSE, gamma = FALSE),
               FLAG_MODEL = list(BASELINE = TRUE, SSE = FALSE, SSI = TRUE,
                                 SSEB = FALSE, SSIB = FALSE),
               FLAG_FILTER = list(end_day = FALSE, tot_infs = TRUE))

#******************
#* 1. SSE
COMP_FOLDER = "~/GitHub/computing/mcmc/sse/"
PLOT_CI_PARAMS(df_sse, COMP_FOLDER, fig_num = '01', PDF = TRUE, 
               FLAG_PARAM = list(r0 = FALSE, k = TRUE, 
                                 alpha = FALSE, gamma = FALSE),
               FLAG_MODEL = list(BASELINE = FALSE, SSE = TRUE, SSI = FALSE,
                                 SSEB = FALSE, SSIB = FALSE),
               FLAG_FILTER = list(end_day = FALSE, tot_infs = TRUE))

#****************
#FIXED
PLOT_CI_PARAMS_FIXED(df_sse1, COMP_FOLDER, fig_num = '001', PDF = TRUE, 
                                 FIXED_PARAM = list(r0 = TRUE, k = FALSE,
                                                    alpha = FALSE, gamma = FALSE),
                                 VAR_PARAM = list(r0 = FALSE, k = TRUE,
                                                  alpha = FALSE, gamma = FALSE),
                                 FLAG_MODEL = list(BASELINE = FALSE, SSE = TRUE, SSI = FALSE,
                                                   SSEB = FALSE, SSIB = FALSE))


#******************
#* 1. SSI 
COMP_FOLDER = "~/GitHub/computing/mcmc/ssi/"
PLOT_CI_PARAMS(df_ssi1, COMP_FOLDER, PDF = FALSE, 
               FLAG_PARAM = list(r0 = FALSE, k = TRUE, 
                                 alpha = FALSE, gamma = FALSE),
               FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = TRUE,
                                 SSEB = FALSE, SSIB = FALSE),
               FLAG_FILTER = list(end_day = FALSE, tot_infs = TRUE))

#SSI FIXED
PLOT_CI_PARAMS_FIXED(df_ssi1, COMP_FOLDER, fig_num = '01', PDF = TRUE, 
                     FIXED_PARAM = list(r0 = TRUE, k = FALSE,
                                        alpha = FALSE, gamma = FALSE),
                     VAR_PARAM = list(r0 = FALSE, k = TRUE,
                                      alpha = FALSE, gamma = FALSE),
                     FLAG_MODEL = list(BASELINE = FALSE, SSE = FALSE, SSI = TRUE,
                                       SSEB = FALSE, SSIB = FALSE))
