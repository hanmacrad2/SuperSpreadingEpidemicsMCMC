#Package install steps

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


