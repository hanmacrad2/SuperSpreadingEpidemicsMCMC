#****************************************************************
# DATA WITH MISSING DATA
#****************************************************************
library(SuperSpreadingEpidemicsMCMC)

DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"

#************************
#SIMULATE DATA
#************************

#BASELINE
file_name = 'data_baseline.rds'
data_ssnb = SIMULATE_EPI_SSNB(num_days = 30)
plot.ts(data_baseline2, main = 'Baseline Simulated data', ylab = 'Daily infections')
saveRDS(data_baseline, paste0(DATA_FOLDER, file_name))

data_baseline = readRDS(paste0(DATA_FOLDER, file_name))
data_baseline_missing_08 = round(0.8*data_baseline)
file_name = 'data_baseline_miss_08.rds'
plot.ts(data_baseline_missing_08, main = 'Baseline Simulated data (20% Missing)', ylab = 'Daily infections')
saveRDS(data_baseline_missing_08, paste0(DATA_FOLDER, file_name))

#RUN MCMC

#DO INFERENCE

#Sub sam