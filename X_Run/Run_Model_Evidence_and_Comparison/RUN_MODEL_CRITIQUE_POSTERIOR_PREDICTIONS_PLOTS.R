#RUN POSTERIOR PREDICTIONS
library(SuperSpreadingEpidemicsMCMC)
ls("package:SuperSpreadingEpidemicsMCMC")

#PROJECT_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/"
#DATA SPECFIC 
DATA = 'NZ_DATA_WAIT_21_SUBSET_I/'
DATA = 'NZ_DATA_WAIT_21/'
DATA = 'BASELINE/'
DATA = 'SSEB/'
DATA = 'SSE/'
OUTER_FOLDER = paste0(PROJECT_FOLDER, DATA)

#DATA

#**********************
#SSE
#**********************
#file_name = 'mcmc_sse_d_15_15.rds'
sim_vals = list(R0 = R0, k = k)
matrix_sim_sse = RUN_POSTERIOR_PREDICTIVE_PLOTS(data_sse, sim_vals,
                                                OUTER_FOLDER, file_name,
                                                SIM_DATA = TRUE,
                                                FLAGS_MODELS = list(BASELINE = FALSE, SSEB = FALSE,
                                                                    SSE = TRUE, SSIB = FALSE, SSI = FALSE))


#1. BASELINE
matrix_sim = RUN_POSTERIOR_PREDICTIVE_PLOTS(data_wait_08_21, OUTER_FOLDER,
                                            SIM_DATA = FALSE,
                               FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE,
                                                   SSE = FALSE, SSIB = FALSE, SSIR = FALSE))

matrix_sim_base = RUN_POSTERIOR_PREDICTIVE_PLOTS(data_baseline, OUTER_FOLDER,
                                            FLAGS_MODELS = list(BASELINE = TRUE, SSEB = FALSE,
                                                                SSE = FALSE, SSIB = FALSE, SSIR = FALSE))

#1. SSEB
matrix_sim_sseb = RUN_POSTERIOR_PREDICTIVE_PLOTS(data_wait_08_21, OUTER_FOLDER,
                                                 SIM_DATA = TRUE,
                                            FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE,
                                                                SSE = FALSE, SSIB = FALSE, SSIR = FALSE))

matrix_sim_sseb = RUN_POSTERIOR_PREDICTIVE_PLOTS(data_wait_08_21, OUTER_FOLDER,
                                                 SIM_DATA = FALSE,
                                                 FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE,
                                                                     SSE = FALSE, SSIB = FALSE, SSIR = FALSE))

matrix_sim_sseb1 = RUN_POSTERIOR_PREDICTIVE_PLOTS(data_wait_08_21_sub1, OUTER_FOLDER,
                                                  SIM_DATA = FALSE,
                                                 FLAGS_MODELS = list(BASELINE = FALSE, SSEB = TRUE,
                                                                     SSE = FALSE, SSIB = FALSE, SSIR = FALSE))
