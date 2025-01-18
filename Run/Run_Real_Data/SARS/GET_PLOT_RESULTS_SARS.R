#***********************
# SARS RESULTS
#***********************

#EPI DATA
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/SARS/'
file_name = 'df_sars_canada_03.rds'
df_sars = readRDS(paste0(RESULTS_FOLDER, file_name))
df_sars$cases = df_sars$total_cases

title = 'SARS Outbreak Canada, 2003. Daily Cases'
data_type = 'sars_canada'
PLOT_EPI_DATA_DATE_PDF(df_sars, RESULTS_FOLDER, title, data_type) 

#*********************************************************
#* GET DATA FILES

RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/SARS/RESULTS/MCMC/'

# List all RDS files in the directory
files <- list.files(path = RESULTS_FOLDER, pattern = "*.rds", full.names = TRUE)

#models
models <- c("baseline", "sse", "ssi", "sseb", "ssib")
model_files <- list()

# Loop through each model to filter files
for (model in models) {
  model_files[[model]] <- files[grepl(model, files)]
}

# Initialize a list to hold the loaded data
data_list <- list()

# Loop through each model and load the files
for (model in models) {
  data_list[[model]] <- lapply(model_files[[model]], readRDS)
}

# Optionally, print the structure of the data_list to verify
str(data_list)

#MCMC
mcmc_baseline = data_list$baseline[[1]]
mcmc_sse = data_list$sse[[1]]
mcmc_ssi = data_list$ssi[[1]]
mcmc_sseb = data_list$sseb[[1]]
mcmc_ssib = data_list$ssib[[1]]

#MODEL COMPARISON
file_name = 'vec_mod_ev2024-05-21_14-18-25_sars_canada_03.rds'
vec_mod_ev = readRDS(paste0(RESULTS_FOLDER, file_name))

file_name = 'post_probs_2024-05-21_14-18-25_sars_canada_03.rds'
post_probs = readRDS(paste0(RESULTS_FOLDER, file_name))

#***********************
#PLOT MCMC RESULTS
RESULTS_FOLDER = '~/GitHub/computing/REAL_DATA/SARS/RESULTS/'

data_title = 'SARS Outbreak Canada, 2003'
main_title = bquote(paste(.(data_title, '. MCMC Results')))
xlimits = c(0, 3.0)

PLOT_MCMC_REAL_DATA(epidemic_data, RESULTS_FOLDER, xlimits,
                    list_mcmc = list(Baseline = mcmc_baseline, SSE = mcmc_sse,
                                     SSI = mcmc_ssi, SSEB = mcmc_sseb, SSIB = mcmc_ssib), MODEL_COLORS,
                    main_title,  plot_margin = c(5.0, 5.2, 4.5, 1.5), cex = 2.45) #bottom left top right
