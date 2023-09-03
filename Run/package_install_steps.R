#Package install steps

#1. STEPS!! 17/07/23
library(devtools)

#1. SET WD to Parent directory
#Code: cd ~/Documents/GitHub

#2. BUILD PACKAGE 
#R CMD build SuperSpreadingEpidemicsMCMC

#3. INSTALL PACKAGE
#Code: R CMD INSTALL SuperSpreadingEpidemicsMCMC_0.1.0.tar.gz


#***************************************************************8
#OPTION 2

# Load the devtools package
library(devtools)

# Build and install the package
outer_folder = "~/GitHub"
setwd(outer_folder)

build("/SuperSpreadingEpidemicsMCMC")
install(SuperSpreadingEpidemicsMCMC)

#LIST OF FUNCTIONS IN PACKAGE
functions_in_package <- ls("package:SuperSpreadingEpidemicsMCMC", pattern = "^[a-zA-Z]", all.names = TRUE)

# Print the list of functions
print(functions_in_package)

#REMOVE OLD VERSION
remove.packages("SuperSpreadingEpidemicsMCMC")












library(SuperSpreadingEpidemicsMCMC)
ls("package:SuperSpreadingEpidemicsMCMC")

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

#install.packages(devtools) *If need to install

#2. library(devtools)

#3. document()

#FINAL: LOAD LIBRARY
# library(brocolors)
# help(plotcolours)

#DEV TOOLS WAY
library(devtools)

#INSTALL PACKAGE
check()
install()

#Rtools (don't need these lines)
#find_rtools() 
#library(devtools)
#install.packages("Rtools")
#library(Rtools)
