#APPLY SSEC_MCMC

#Plot new
R0X = 1.2; kX = 0.16
plot.new()
par(mfrow = c(3, 3))

#Simulated data
ssec_data = SIMULATE_EPI_SSEC()
#Plot
plot.ts(ssec_data,
        main = paste0('Simulated SSEC data, R0: ', R0X, ' k: ', kX),
        ylab = 'Daily infection count')

#2. Simulated data: RANGE OF R0
range_r0 = seq(from = 0.7, to = 2.2, by = 0.1)
range_r0

#PLOT
plot.new()
par(mfrow = c(4, 4))

#LOOP VALUES
for (i in seq_along(range_r0)){
  print(i)
  #Data
  ssec_data = SIMULATE_EPI_SSEC(R0 = range_r0[i])
  #Plot
  plot.ts(ssec_data,
          main = paste0('Simulated SSEC data, R0: ', range_r0[i], ' k: ', kX),
          ylab = 'Daily infection count')
}

#2. Simulated data: RANGE OF k
range_k1 = c(rep(0.01, 4), rep(0.05, 4), rep(0.1, 4), rep(0.25, 4))
range_k1

range_k2 = c(rep(0.5, 4), rep(0.75, 4), rep(1.0, 4), rep(5, 4))
range_k2

range_k2 = rep(c(0.5, 0.75, 1.0, 5), 4)
range_k2

#PLOT
plot.new()
par(mfrow = c(4, 4))

#LOOP VALUES
for (i in seq_along(range_k2)){
  print(i)
  #Data
  ssec_data = SIMULATE_EPI_SSEC(num_days = 110, k = range_k2[i])
  #Plot
  if (i == 1){
    plot.ts(ssec_data,
            main = bquote(bold("SSEC data ~ NegBin(R0 * " ~ lambda ~ ", k). " ~
            k: ~ .(range_k2[i]) ~ ". R0: " ~ .(R0X))),
            ylab = 'Daily infection count')
  } else {
    plot.ts(ssec_data,
            main = bquote(bold(k: ~ .(range_k2[i]) ~ ". " ~ R0: ~ .(R0X))),
            ylab = 'Daily infection count')
  }
}
