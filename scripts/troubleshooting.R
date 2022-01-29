# make demographic distributions for species

library (RCurl)
getURL()

rm(list=ls())
source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')

source('scripts/PW_functions_local.R')

# Load the per-year data (without aggregating branches)
indv.data.2013 <- get.indv.data(year = 2013,branches=T)
#head(indv.data.2013)

indv.data.2018 <- get.indv.data(year = 2018,branches=T)
#head(indv.data.2018)

indv.data.2019 <- get.indv.data(year = 2019,branches=T)
#head(indv.data.2019)

indv.data.2020 <- get.indv.data(year = 2020,branches=T)
head(indv.data)

str(indv.data.2020)

!(1000:5153 %in% indv.data.2013$Num)

# THe following lines identify which plots dupolicated numbers are in
indv.data <- indv.data.2018

# Check to see if tag numbers are duplicated in a particular year
length(unique(indv.data$Num)) # How many UNIQUE tag numbers
length(indv.data$Num) # How many TOTAL tag numbers
# Make tables with the duplicated number and which plot it's in
first <- data.frame(onetag=indv.data$Num[duplicated(indv.data$Num,fromLast=TRUE)],oneplot=indv.data$Plot[duplicated(indv.data$Num,fromLast=TRUE)])
second <- data.frame(twotag=indv.data$Num[duplicated(indv.data$Num)],twoplot=indv.data$Plot[duplicated(indv.data$Num)])
# Sort the tables so that the numbers are in the same order
first <- first[order(first$onetag),]
second <- second[order(second$twotag),]
# Bind the first and sxwecond dup tables
dups <- cbind(first,second)
# Add a check to see if duplicates are in the same plot
dups$sameplot <- apply(dups, 1, FUN=function(x) x[2] == x[4])
dups[order(dups$sameplot),]


grep("DUP",indv.data.2018$Notes)
indv.data.2018[1559,]
indv.data.2018[37,]



#checking which plot a particular tag number is in - fixed this particular error in CSV
indv.data.2013$Plot[indv.data.2013$Num==9196.3]

# Identify new tag numbers in 2018 that weren't in the 2013 data
# Including ones less than 5154 that must have been out
newtags <- sort(indv.data.2018$Num[!(indv.data.2018$Num %in% indv.data.2013$Num)])
newtags <- newtags[newtags < 5500]

indv.data.2018[indv.data.2018$Num %in% newtags,]
