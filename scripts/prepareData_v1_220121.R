# make demographic distributions for species
rm(list=ls())
update.packages(c('Rcurl','data.table','ape','picante','vegan','permute'))
library(RCurl)

source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')

source('scripts/PW_functions_local.R')

# Load the per-year data (without aggregating branches)
indv.data.2013 <- get.indv.data(year = 2013,branches=F)
dim(indv.data.2013)
head(indv.data.2013)

#indv.data.2013.B <- get.indv.data(year = 2013,branches=T)
#dim(indv.data.2013.B)
#head(indv.data.2013.B)

indv.data.2018 <- get.indv.data(year = 2018,branches=F)
dim(indv.data.2018)
head(indv.data.2018)

#indv.data.2018.B <- get.indv.data(year = 2018,branches=T)
#dim(indv.data.2018.B)
#head(indv.data.2018.B)

indv.data.2019 <- get.indv.data(year = 2019,branches=F)
dim(indv.data.2019)
head(indv.data.2019)

#indv.data.2019.B <- get.indv.data(year = 2019,branches=T)
#dim(indv.data.2019.B)
#head(indv.data.2019.B)

# This one gives an error - skip for now

plot.list.2020 <- get.plot(2020)
# now remove plots not sampled
plot.list.2020 <- plot.list.2020[-c(5,6,20,25,28,30,40:42,46,49,51:54)]
plot.list.2020
indv.data.2020 <- get.indv.data(year = 2020,branches=F,plot.list=plot.list.2020)
dim(indv.data.2020)
head(indv.data.2020)
# str(indv.data.2020)

## which numbers were missing in 2013
head(sort(indv.data.2013$Num))
tail(sort(indv.data.2013$Num))
maxnum.2013 <- max(indv.data.2013$Num)

# There's one individual missing a number! Where is that?
indv.data.2013[which(indv.data.2013$Num==-999),]

# Setting that one aside, which numbers are missing between 1001 and 5153 in 2013
allNum13 <- 1001:5153
missNum13 <- which(!(allNum13 %in% indv.data.2013$Num))
allNum13[missNum13]
length(allNum13[missNum13]) # number of missing tag numbers

# Are there duplicated tags in any year?
# put the three individual data dataframes into a list so that we can loop through them and do operations on them
all.indv.data <- list(indv.data.2013,indv.data.2018,indv.data.2019,indv.data.2020)
for (i in 1:length(all.indv.data)) 
{
  indv.data <- indv.data.2013
  tagT <- table(indv.data.2013$Num)
  print(c(i,tagT[which(tagT>1)]))
}
## APPARENTLY NO DUPLICATED TAGS IN ANY YEAR

# Identify new tag numbers in 2018 that weren't in the 2013 data

# Including ones less than 5154 that must have been out
newtags <- sort(indv.data.2018$Num[!(indv.data.2018$Num %in% indv.data.2013$Num)])
newtags
(newtags.lt.maxnum <- newtags[newtags <= maxnum.2013])
(newtags.gt.maxnum <- newtags[newtags > maxnum.2013])
maxnum.2018 <- 6325 # entered manually due to two outliers

# those last two don't look right - worth fixing in csv
indv.data.2018[which(indv.data.2018$Num==27775),]
indv.data.2018[which(indv.data.2018$Num==40070),]

# Identify new tag numbers in 2019 that weren't in the 2018 data
# Including ones less than 5154 that must have been out
newtags.19 <- sort(indv.data.2019$Num[!(indv.data.2019$Num %in% indv.data.2018$Num)])
newtags.19
(newtags.19[newtags.19 <= maxnum.2018])
(newtags.19[newtags.19 > maxnum.2018])
maxnum.19 <- 6951

# those last five don't look right - worth fixing in csv
indv.data.2019[which(indv.data.2019$Num==9595),]
indv.data.2019[which(indv.data.2019$Num==20804),]
indv.data.2019[which(indv.data.2019$Num==21223),]
indv.data.2019[which(indv.data.2019$Num==47990),]
indv.data.2019[which(indv.data.2019$Num==48070),]

# Identify new tag numbers in 2020 that weren't in the 2019 data
# Including ones less than 5154 that must have been out
newtags.20 <- sort(indv.data.2020$Num[!(indv.data.2020$Num %in% indv.data.2019$Num)])
newtags.20
(newtags.20[newtags.20 <= maxnum.19])
(newtags.20[newtags.20 > maxnum.19])

# there's one funny one
indv.data.2020[which(indv.data.2020$Num==282873),]

## Now make a list of all numbers that appear across all years
# The loop here makes this so it will run smoothly if additional years are added
allNums <- c()

# all.indv.data is the list made above, where each item is one years individual data. How many years does it have:
length(all.indv.data)

# make an empty variable, and then loop through the individual data files and append all the numbers end to end
allNums <- c()
for (i in 1:length(all.indv.data)) allNums <- c(allNums,all.indv.data[[i]]$Num)
head(allNums)
# now reduce to the unique ones
allNums <- sort(unique(allNums))
length(allNums)
head(allNums)
tail(allNums)

##### edit from here for 2020

# there are some numbers above that should be fixed. In the meantime, let's make a dataframe with all individuals across all years, which we can use to start checking problems
allIndv <- data.frame(Num=allNums,P13=NA,P18=NA,P19=NA,P20=NA,S13=NA,S18=NA,S19=NA,S20=NA)

# now use the match command to match up the plot for each number in each year, and assign it to the right row
Pn <- c('P13','P18','P19','P20')
Sn <- c('S13','S18','S19','S20')
i=1
for (i in 1:length(all.indv.data))
{
  y2a <- match(allIndv$Num,all.indv.data[[i]]$Num)
  allIndv[,Pn[i]] <- all.indv.data[[i]]$Plot[y2a]
  allIndv[,Sn[i]] <- all.indv.data[[i]]$Species[y2a]
}
head(allIndv)

## This data.frame may be quite useful for quickly identifying problem individuals. For example, what numbers were present in plot 1345 in 2018 that were missed in 2013. Maybe that's the right number for the -999 individual
allIndv[which(allIndv$P18=='PPW1345'),c('Num',Pn)]
# Indv 2276 is missing in 2013 and present in later years. Is that the right number for the one with -999?

# Now, let's see which numbers were assigned in different plots in different years. To do this we'll use the apply function which can apply a function either across the rows or columns of a matrix or dataframe. We'll make a new function which counts the number of unique entries. So if all three plot IDs were the same for a number, the function returns 1. If it's more than 1, indicates a problem
lunique <- function(x) 
{
  nonNAx <- x[which(!is.na(x))]
  return(length(unique(nonNAx)))
}

allIndv$nP <- apply(allIndv[,Pn],1,lunique)
head(allIndv$nP)
(multPlots <- which(allIndv$nP>1))

# OK, 25 individuals that moved between plots! We need to fix those. Here they are:
allIndv[multPlots,]

# and which ones change species
allIndv$nSp <- apply(allIndv[,Sn],1,lunique)
head(allIndv$nSp)
(multSp <- which(allIndv$nSp>1))

# OK, almost 300 individuals with more than one Sp ID! We'll need to fix or exclude these. Here they are:
allIndv[multSp,c('Num',Sn)]

## which individuals are missing from 2018 and present in 2013 and 2019
midNA <- function(x)
{
  if (!is.na(x[1]) & is.na(x[2]) & !is.na(x[3])) ret <- 1 else ret <- 0
  return(ret)
}
allIndv[apply(allIndv[,Pn[1:3]],1,midNA)==1,]
allIndv[apply(allIndv[,Pn[2:4]],1,midNA)==1,]

########
# OLD CODE BELOW HERE
# 
# # THe following lines identify which plots dupolicated numbers are in
# 
# # Check to see if tag numbers are duplicated in a particular year
# length(unique(indv.data$Num)) # How many UNIQUE tag numbers
# length(indv.data$Num) # How many TOTAL tag numbers
# # Make tables with the duplicated number and which plot it's in
# first <- data.frame(onetag=indv.data$Num[duplicated(indv.data$Num,fromLast=TRUE)],oneplot=indv.data$Plot[duplicated(indv.data$Num,fromLast=TRUE)])
# second <- data.frame(twotag=indv.data$Num[duplicated(indv.data$Num)],twoplot=indv.data$Plot[duplicated(indv.data$Num)])
# # Sort the tables so that the numbers are in the same order
# first <- first[order(first$onetag),]
# second <- second[order(second$twotag),]
# # Bind the first and sxwecond dup tables
# dups <- cbind(first,second)
# # Add a check to see if duplicates are in the same plot
# dups$sameplot <- apply(dups, 1, FUN=function(x) x[2] == x[4])
# dups[order(dups$sameplot),]
# 
# 
# grep("DUP",indv.data.2018$Notes)
# indv.data.2018[1559,]
# indv.data.2018[37,]
# 
