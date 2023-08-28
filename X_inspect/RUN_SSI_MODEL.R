#SSI MODEL INSPECT
par(mfrow=c(4,4))
num_days = 50

#k = 0.2, 0.5, 0.9, 1.2
#ssi_matrix = matrix(0, nrow = num_days, ncol = 16)
R0 = 1.5
k = 1.0 #0.1 0.3 #9 #1.5 
data_ssi = SIMULATE_EPI_SSI(num_days = num_days, R0 = R0, k = k)
plot.ts(data_ssi$epidemic_data,
        main = paste0('SSI Data, R0: ', R0, ', k: ', k),
        ylab = 'Daily count')