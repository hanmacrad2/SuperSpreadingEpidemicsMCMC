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

#COUNTIES MANAUKA
file_name = 'data_cm_tot_aug_21.csv'

data_file_cm = read.csv(paste0(DATA_FOLDER, file_name))
data_cm_08_21 = data_file_cm$Cases
plot.ts(data_cm_08_21, ylab = 'Infection count', main = ' Counties Manauka, August 2021')

#SUBSET ONE
data_cm_08_21_sub1 = data_cm_08_21[3:22]
plot.ts(data_cm_08_21_sub1, ylab = 'Infection count', main = ' Counties Manauka, August 2021')
