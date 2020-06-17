# the sp library has all the interpolation stuff in R
library(sp)

# set the working directory in R Studio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load the stations
stations <- read.table("niederschlagsdaten/nieder_stationen.txt", 
                       header = TRUE, sep = ";", dec = ".", fileEncoding="latin1")
stations #check

# drop the columns we don't need:
stations <- stations[c('Name', 'east', 'north')]
stations #check

# load the daily precipitation data
precip <- read.table("niederschlagsdaten/nieder_tageswerte.txt", 
                     header = TRUE, sep = ";", dec = ".", fileEncoding="latin1")
precip # you get the idea...

# we will only look at June 14, 2020 - that day has a nice distribution of different 
# values across NRW
precip <- precip[(precip$Datum == "2020-06-14"),]
precip 

# only keep the value and the station name
precip <- precip[c('Name', 'Tagessumme')] 
precip

# join on the station name to get the coordinates
precip <- merge(stations, precip, by="Name")
precip

# write that out, so we can look it in a GIS etc.
write.csv(precip, "June14.csv") 


