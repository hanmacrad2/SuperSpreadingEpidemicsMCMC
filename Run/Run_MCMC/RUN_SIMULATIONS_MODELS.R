#DATA FROM MODEL SIMULATIONS
library(SuperSpreadingEpidemicsMCMC)
ls("package:SuperSpreadingEpidemicsMCMC")

#PARAMETERS
R0X = 1.6
par(mfrow = c(2,2))

#1. BASELINE
model_type = 'Baseline'
data_base3 = SIMULATE_EPI_BASELINE(R0X)
plot.ts(data_base, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' simulated data. R0 = ', R0X))

#2. SSNB
model_type = 'SS-NB'
data_ssnb3 = SIMULATE_EPI_SSNB(R0X)
plot.ts(data_base, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated data. R0: ', R0X))

#3. SSEB
model_type = 'SSEB'
data_sseb3 = SIMULATE_EPI_SSEB(alphaX = 0.6, betaX = 0.1) #*Should print want R0x is 
plot.ts(data_sseb, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated Data. R0: ', R0X))

#4. SSIB
model_type = 'SSEB'
data_sseb3 = SIMULATE_EPI_SSIB(alphaX = 0.6, betaX = 0.1) #*Should print want R0x is 
plot.ts(data_sseb, ylab = 'Infection count', xlab = 'Day',
        main = paste0(model_type, ' Simulated Data. R0: ', R0X))

#5.
model_type = 'SS_IR'
data_ssir = SIMULATE_EPI_SSIR(1.6)
