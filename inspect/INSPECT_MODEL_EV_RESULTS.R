#INSPECT MODEL EVIDENCE RESULTS
#*MAKE FUNCTION

#LOG of Bayes Factors
log_bfs1 = 0.5*ests_base - 0.5*ests_base - 0.25*ests_sseb - 0.25*ests_ssib

#Log bfs (recipricol)
log_bfs2 = 0.25*ests_sseb + 0.25*ests_ssib - 0.5*ests_base- 0.25*ests_sseb - 0.25*ests_ssib

#LOG PLOTS
plot(seq_along(log_bfs2), log_bfs2,
     ylim = c(min(log_bfs2)-0.05, max(log_bfs2)+0.05), #lwd = 1, pch = 19,
     main = 'Ratio of Model evidences via Importance sampling. Base data. Equal priors NS vs SS',
     xlab = 'rep', ylab = 'Ratio of Model evidences')

#BAYES FACTORS (50:50 weight on regular spread vs superspread)
bfs1 = (0.5*ests_base)/(0.5*ests_base + 0.25*ests_sseb + 0.25*ests_ssib)

#PLOTS
plot(seq_along(bfs1), bfs1,
     ylim = c(min(bfs1)-0.05, max(bfs1)+0.05), #lwd = 1, pch = 19,
     main = 'Ratio of Model evidences via Importance sampling. Base data. Equal priors NS vs SS',
     xlab = 'rep', ylab = 'Ratio of Model evidences')

#PLOTS
plot(seq_along(bfs1), bfs1,
     ylim = c(min(bfs1)-0.05, max(bfs1)+0.05), #lwd = 1, pch = 19,
     main = 'Ratio of Model evidences via Importance sampling. Base data. Equal priors NS vs SS',
     xlab = 'rep', ylab = 'Ratio of Model evidences')

#********************
#PLOTS (RECIPROCAL)
bfs2 = 1/bfs1 #(0.5*ests_base)/(0.5*ests_base + 0.25*ests_sseb + 0.25*ests_ssib)
plot(seq_along(bfs2), bfs2,
     ylim = c(min(bfs1)-0.05, max(bfs1)+0.05), #lwd = 1, pch = 19,
     main = '(Reciprocal)Ratio of Model evidences via Importance sampling. Base data. Equal priors NS vs SS',
     xlab = 'rep', ylab = 'Ratio of Model evidences')

#************************
#* MODEL EVIDENCE 
#PLOTS
par(mfrow = c(2,1))
plot(seq_along(ests_base), ests_base,
     ylim = c(min(ests_base)-50, max(ests_base)+50), lwd = 1, pch = 19,
     main = 'Model evidence (Baseline) via Importance sampling. Base data',
     xlab = 'rep', ylab = 'Model evidence estimate')

#PLOT 2 (Close up)
plot(seq_along(ests_base), ests_base,
     ylim = c(min(ests_base)-0.5, max(ests_base)+0.5), lwd = 1, pch = 19,
     main = 'Model evidence (Baseline) via Importance sampling. Base data. (Close up)',
     xlab = 'rep', ylab = 'Model evidence estimate')

#************
#* SSEB & SSIB
#************

#PLOTS
par(mfrow = c(2,1))
plot(seq_along(ests_ssib), ests_ssib,
     ylim = c(min(ests_ssib)-50, max(ests_ssib)+50), lwd = 1, pch = 19,
     main = 'Model evidence (SSIB) via Importance sampling. Base data',
     xlab = 'rep', ylab = 'Model evidence estimate')

#PLOT 2 (Close up)
plot(seq_along(ests_ssib), ests_ssib,
     ylim = c(min(ests_ssib)-0.5, max(ests_ssib)+0.5), lwd = 1, pch = 19,
     main = 'Model evidence (SSIB) via Importance sampling. Base data. (Close up)',
     xlab = 'rep', ylab = 'Model evidence estimate')