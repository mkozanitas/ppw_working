# RUN prepareData again if csv's have been
rm(list=ls())

all.id <- readRDS('data/allid.Rdata')
all.idb <- readRDS('data/allidb.Rdata')

# id13 <- all.id[[1]]
# id18 <- all.id[[2]]
# id19 <- all.id[[3]]
# id20 <- all.id[[4]]
# head(id13)

years <- c(2013,2018,2019,2020)
# etc.

# Print dups for each year. And for now, remove dups before checking for other problems, e.g. moving between plots, species, etc.
i=1
for (i in 1:4)
{
  print(years[i])
  nNum <- table(all.id[[i]]$Num)
  print(nNum[which(nNum>1)])
  dups <- as.numeric(names(nNum[which(nNum>1)]))
  if (length(dups)>0) all.id[[i]] <- all.id[[i]][-which(all.id[[i]]$Num %in% dups)]
}

for (i in 1:4) print(tail(sort(all.id[[i]]$Num)))

### ONE BAD TAG: 39332 in 2018
all.id[[2]][all.id[[2]]$Num==39332,]

# all.id is a list made above, where each item is one years individual data. How many years does it have:
length(all.id)

# make an empty variable, and then loop through the individual data files and append all the numbers end to end
allNums <- c()
for (i in 1:length(all.id)) allNums <- c(allNums,all.id[[i]]$Num)
head(sort(allNums))
tail(sort(allNums))

# now reduce to the unique ones
allNums <- sort(unique(allNums))
length(allNums)
head(allNums)
tail(allNums)

# there are some numbers above that should be fixed. In the meantime, let's make a dataframe with all individuals across all years, which we can use to start checking problems
allIndv <- data.frame(Num=allNums,P13=NA,P18=NA,P19=NA,P20=NA,S13=NA,S18=NA,S19=NA,S20=NA)

# now use the match command to match up the plot for each number in each year, and assign it to the right row
Pn <- c('P13','P18','P19','P20')
Sn <- c('S13','S18','S19','S20')
i=1
for (i in 1:length(all.id))
{
  y2a <- match(allIndv$Num,all.id[[i]]$Num)
  allIndv[,Pn[i]] <- all.id[[i]]$Plot[y2a]
  allIndv[,Sn[i]] <- all.id[[i]]$Species[y2a]
}
head(allIndv)

## This data.frame may be quite useful for quickly identifying problem individuals. 

# Now, let's see which numbers were assigned in different plots in different years. To do this we'll use the apply function which can apply a function either across the rows or columns of a matrix or dataframe. We'll make a new function which counts the number of unique entries. So if all three plot IDs were the same for a number, the function returns 1. If it's more than 1, indicates a problem
lunique <- function(x) 
{
  nonNAx <- x[which(!is.na(x))]
  return(length(unique(nonNAx)))
}

allIndv$nP <- apply(allIndv[,Pn],1,lunique)
head(allIndv$nP)
(multPlots <- which(allIndv$nP>1))

# OK, 19 individuals that moved between plots (not counting ones that were duplicated in one or more years, which would need to resolved first) We need to fix those. Here they are:
allIndv[multPlots,1:9]

# and which ones change species
allIndv$nSp <- apply(allIndv[,Sn],1,lunique)
head(allIndv$nSp)
(multSp <- which(allIndv$nSp>1))

# what did we find?
table(allIndv$nSp)



# Oh, interesting - individuals not assigned to any species? And some assigned to different species every year? Let's look at 0, 3
probs <- which(allIndv$nSp %in% c(0,3,4))
allIndv[probs,]

# OK, almost 300 individuals with more than one Sp ID! We'll need to fix or exclude these. Here they are:
head(allIndv[multSp,1:9])

## how many are the same individuals
table(allIndv$nP,allIndv$nSp)

#make combined list
allMults <- allIndv[union(multPlots,multSp),]
write.csv(allMults[order(allMults$nSp,allMults$nP,decreasing = T),],'data/mult-plots-species.csv')

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