#NZ DATA 

#*************************
# NZ OUTBREAK AUGUST 2021
DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
file_name =  'data_nz_tot_aug_21.csv'

data_file_nz = read.csv(paste0(DATA_FOLDER, file_name))
data_nz_08_21 = data_file_nz$Cases
plot.ts(data_nz_08_21, ylab = 'Infection count', main = ' NZ, August 2021')

#SUBSET ONE
data_nz_08_21_sub1 = data_nz_08_21[1:22]
plot.ts(data_nz_08_21_sub1, ylab = 'Infection count', main = ' NZ, August 2021')

#****************
#WAITMETA
#****************

DATA_FOLDER = "~/GitHub/SuperSpreadingEpidemicsMCMC/data/"
file_name =  'data_waitemata_aug_21.csv'

data_file_waitemata = read.csv(paste0(DATA_FOLDER, file_name))
data_wait_08_21 = data_file_waitemata$Cases
data_wait_08_21 = data_wait_08_21[2:length(data_wait_08_21)]
plot.ts(data_wait_08_21, ylab = 'Infection count', main = 'Waitemata, August 2021')

#SUBSET ONE
t1 = 2; t2 = 11
data_wait_08_21_sub1 = data_wait_08_21[t1:t2]
plot.ts(data_wait_08_21_sub1, ylab = 'Infection count',
        main = paste0('Waitemata, August 2021. Days ', t1, ' to ', t2, ' (Subset 1)'))

#SUBSET TWO
t1 = 10; t2 = 20
data_wait_08_21_sub2 = data_wait_08_21[t1:t2]
plot(c(10:20), type = 'l',
  data_wait_08_21_sub2, xlab = 'Time', ylab = 'Infection count',
        main = paste0('Waitemata, August 2021. Days ', t1, ' to ', t2, ' (Subset II)'))


#COUNTIES MANAUKA
file_name = 'data_cm_tot_aug_21.csv'

data_file_cm = read.csv(paste0(DATA_FOLDER, file_name))
data_cm_08_21 = data_file_cm$Cases
plot.ts(data_cm_08_21, ylab = 'Infection count', main = ' Counties Manauka, August 2021')

#SUBSET ONE
data_cm_08_21_sub1 = data_cm_08_21[3:22]
plot.ts(data_cm_08_21_sub1, ylab = 'Infection count', main = ' Counties Manauka, August 2021')
