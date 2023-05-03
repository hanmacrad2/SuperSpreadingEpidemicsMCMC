#DATA FROM MODEL SIMULATIONS
library(SuperSpreadingEpidemicsMCMC)
ls("package:SuperSpreadingEpidemicsMCMC")
DATA_FOLDER = '~/GitHub/SuperSpreadingEpidemicsMCMC/data/'

#PARAMS
num_models = 5; num_days = 50
matrix_data = matrix(NA, num_models, num_days)

#PLOT
plot_sim_model <- function(epi_data, model_type, R0X){
  
  #PLOT
  plot(seq_along(epi_data), epi_data, type = 'l',
          ylab = 'Infection count', xlab = 'Day',
          mgp = c(3, 1, 0),
          cex.lab = 1.7, cex.axis = 1.4, cex.main=1.5, #cex.sub= 1.5,
          main = paste0(model_type, ' simulated data. R0 = ', R0X))
  
}
#PLOT ALL DATA
plot_sim_matrix_dat <- function(matrix_data, R0X = 1.6){
  
  #PLOT
  par(mfrow = c(2,3))
  model_type = c('Baseline', 'SSE-NB', 'SSE-B', 'SSI-B', 'SS-IR')
  num_datasets = dim(matrix_data)[1]
  
  for (i in c(1:num_datasets)){
    print(i)
    plot_sim_model(matrix_data[i,], model_type[i],  R0X = R0X)
  }
}

#LOOP & PLOT
plot_sim_matrix_dat(matrix_data)

#PARAMETERS
R0X = 1.6
par(mfrow = c(2,3))

#1. BASELINE
model_type = 'Baseline'
#data_base3 = SIMULATE_EPI_BASELINE(R0X)
matrix_data[1,] = data_base
plot.ts(data_base, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' simulated data. R0 = ', R0X))

#2. SSNB
model_type = 'SS-NB'
data_ssnb = SIMULATE_EPI_SSNB(R0X)
plot.ts(data_ssnb, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated data. R0: ', R0X))
matrix_data[2,] = data_ssnb

#3. SSEB
model_type = 'SSEB'
#data_sseb3 = SIMULATE_EPI_SSEB(alphaX = 0.6, betaX = 0.1) #*Should print want R0x is 
plot.ts(data_sseb, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated Data. R0: ', R0X))
matrix_data[3,] = data_sseb

#4. SSIB
model_type = 'SSIB'
#data_ssib = SIMULATE_EPI_SSIB() #*Should print want R0x is 
plot.ts(data_ssib, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated Data. R0: ', R0X))
matrix_data[4,] = data_ssib

#5.SS-IR or SS-Ind
model_type = 'SS_IR'
data_ssir4 = SIMULATE_EPI_SSIR(1.6)
plot.ts(data_ssir4$epidemic_data, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated Data. R0: ', R0X))
matrix_data[5,] = data_ssir4$epidemic_data

#SAVE THE DATA
saveRDS(matrix_data, file = paste0(DATA_FOLDER, 'matrix_simulated_datasets.rds'))
