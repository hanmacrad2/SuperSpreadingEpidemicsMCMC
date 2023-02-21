#Contact tracing Data
#data_link = "https://github.com/dcadam/covid-19-sse/blob/master/data/secondary_cases.csv"

#Data
data_set = 'secondary_cases.csv'
data_sec_cases <- read.csv(data_set, header = TRUE, sep = ",", dec = ".")

#Plot
hist(data_sec_cases[["secondary.cases"]], 
     breaks = 50,
     main = 'Contact tracing data (secondary cases) Hong Kong',
     xlab = 'Number of secondary cases')

#Mean
#Mean = 1.69697

