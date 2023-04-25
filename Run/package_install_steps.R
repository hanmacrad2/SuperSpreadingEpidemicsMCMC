#Package install steps

library(devtools)
install.packages("Rtools")
library(Rtools)

#1. SET WD to Parent directory
#Code: cd ~/Documents/GitHub

#2. BUILD PACKAGE 
#Code: R CMD build SuperSpreadingEpidemicsMCMC

#3. INSTALL PACKAGE
#Code: R CMD INSTALL SuperSpreadingEpidemicsMCMC_0.1.0.tar.gz

library(SuperSpreadingEpidemicsMCMC)

#NB!!!
#Note: Need to delete unused functions in man and delete Namespace and re-build


#PART 1: TERMINAL INSTALL
#1. Set wd to parent directory

#2. R CMD build brocolors

#3. R CMD INSTALL brocolors_0.1.tar.gz


#PART 2:R Console

# library(brocolors)

#Run functions;
# brocolors() #Didnt work for me without rendering documentation
# plot_crayons()

#PART 3: Document

#1. setwd(parent_directory)

#* install.packages(devtools) *If need to install

#2. library(devtools)

#3. document()

#FINAL: LOAD LIBRARY
# library(brocolors)
# help(plotcolours)


#R TRICKS

#List of named vectors
#USE = sign, not ->
i = 4000
thinning_factor = 10; burn_in_start = 0.2*n_mcmc
if (i%%thinning_factor == 0 & i >= burn_in_start) {
  print('y')
}
