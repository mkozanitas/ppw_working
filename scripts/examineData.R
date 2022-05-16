# RUN prepareData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R
# data files are not synced by git, so prepareData does need to be run locally
rm(list=ls())

# read in list of 4 items, each item with the full data file for all trees in a given year
years <- c(2013,2018,2019,2020)
all.id <- readRDS('data/allid.Rdata')

# these are data files with 'points' (branches) output individually. Not currently being used, and the lines to create these files may be commented out in prepareData
#all.idb <- readRDS('data/allidb.Rdata')


# Print dups for each year. And for now, remove dups before checking for other problems, e.g. moving between plots, species, etc.
dups <- list()
i=1
for (i in 1:4)
{
  print(years[i])
  nNum <- table(all.id[[i]]$Num)
  print(nNum[which(nNum>1)])
  #next line gives you tree numbers that are duplicated
  dN <- as.numeric(names(nNum[which(nNum>1)]))
  dups[[i]] <- all.id[[i]][which(all.id[[i]]$Num %in% dN),]
  dups[[i]] <- dups[[i]][order(dups[[i]]$Num),]
  
  #the next line removes duplicated tags
  if (length(dN)>0) all.id[[i]] <- all.id[[i]][-which(all.id[[i]]$Num %in% dN),]
}

dups[[1]]
dups <- do.call(rbind,dups)
dim(dups)
head(dups)
write.csv(dups,'data/duplicates.csv')

# check that the maximum numbers from each year of survey don't have typos or bad values
for (i in 1:4) print(tail(sort(all.id[[i]]$Num)))

# Pull out all individuals with NA for bSprout - this should always be filled out. As of 5/15/22, no problems identified at this step!!
bNA <- all.id[[1]][which(is.na(all.id[[1]]$bSprout)),]
for (i in 2:4) bNA <- rbind(bNA,all.id[[i]][which(is.na(all.id[[i]]$bSprout)),])
dim(bNA)
head(bNA)
write.csv(bNA,'data/bSprout-NAs.csv')

## Seems there are problems with gCrown - checking here
summary(all.id[[1]]$gCrown) # GOOD
summary(all.id[[2]]$gCrown) # GOOD
summary(all.id[[3]]$gCrown) # WHO is -Inf?
rsel <- which(all.id[[3]]$gCrown==-Inf)
all.id[[3]][rsel,]

summary(all.id[[4]]$gCrown) # Character! Let's make it numeric
all.id[[4]]$gCrown <- as.numeric(all.id[[4]]$gCrown)
summary(all.id[[4]]$gCrown) # Good!

# In post-fire years, check individuals scored as any combination of DEAD & TOPKILL, DEAD & GREEN, TOPKILL & GREEN

# catenate values to see patterns
catVals <- function(x) {
  res <- c()
  for (i in 1:length(x)) res <- paste(res,x[i],sep='')
  return(res)
}

## ONLY DO THIS FOR 2018 AND BEYOND (i in 2:4)
# These are patterns of values for 8 fields in the data (see below) which represent the 'legal' combinations. Any tree that doesn't follow this pattern suggests either a data ehtry or a coding error requiring further investigation. As of 5/15/22 there are no problem saplings, and just 7 and 24 problem trees in 2018 and 2019. But there are >1000 in 2020, so there's some deeper problem we need to figure out.
SA.patts <- c('00NANA1010','01NANA0110','10NANA0101','11NANA0101')
TR.patts <- c('00001010','01000110','10010101','10100101','10110101','11010101','11100101','11110101')

i=4
for (i in 2:4) {
  all.id[[i]]$pattern <- apply(all.id[[i]][,c("Survival","bSprout","Epicormic","Apical","Dead","Live","Topkill","gCrown")],1,catVals)
  print(table(all.id[[i]]$pattern[all.id[[i]]$Type=='SA']))
  badSAs <- which(all.id[[i]]$Type=='SA' & !all.id[[i]]$pattern %in% SA.patts)
  length(badSAs)
  
  print(table(all.id[[i]]$pattern[all.id[[i]]$Type=='TR']))
  badTRs <- which(all.id[[i]]$Type=='TR' & !all.id[[i]]$pattern %in% TR.patts)
  length(badTRs)
  print(all.id[[i]][badTRs,c('Plot','Num','pattern')])
}

# examine trees with particular problem patterns
i=4
all.id[[i]][which(all.id[[i]]$pattern=='11NANA010NA'),]

# more troubleshooting code - commented out for now
# SA18 <- all.id[[2]][which(all.id[[2]]$Type=='SA'),]
# table(TR18$Live,TR18$Topkill)
# table(TR18$Live,TR18$gCrown)
# table(TR18$Topkill,TR18$gCrown)
# 
# table(SA18$Live,SA18$Topkill)
# table(SA18$Live,SA18$gCrown)
# table(SA18$Topkill,SA18$gCrown)

# Now create some combined states, again to look for 'illegal' data combinations. Need to revisit this - some of the '2s' may suggest problems, but not sure.
i=4
for (i in 2:4) {
  all.id[[i]]$Dead <- 1 - all.id[[i]]$Live
  all.id[[i]]$DT <- all.id[[i]]$Dead + all.id[[i]]$Topkill
  all.id[[i]]$TG <- all.id[[i]]$Topkill + all.id[[i]]$gCrown
  all.id[[i]]$DG <- all.id[[i]]$Dead + all.id[[i]]$gCrown
  print(c(years[i],'Dead'))
  print(table(all.id[[i]]$Dead))
  print(c(''))
  print(c(years[i],'Topkill'))
  print(table(all.id[[i]]$Topkill))
  print(c(''))
  print(c(years[i],'gCrown'))
  print(table(all.id[[i]]$gCrown))
  print(c(years[i],'DT'))
  print(table(all.id[[i]]$DT))
  print(c(years[i],'TG'))
  print(table(all.id[[i]]$TG))
  print(c(years[i],'DG'))
  print(table(all.id[[i]]$DG))
}


# all.id is a list made above, where each item is one years individual data. How many years does it have:
length(all.id)

# make an empty variable, and then loop through the individual data files and append all the numbers end to end
allNums <- c()
for (i in 1:length(all.id)) allNums <- c(allNums,all.id[[i]]$Num)

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

# Now check all names across all years
head(allIndv[,c('S13','S18','S19','S20')])
allNames <- sort(c(allIndv$S13,allIndv$S18,allIndv$S19,allIndv$S20))
table(allNames)

write.csv(sort(unique(allNames)),'data/all-spp-names.csv')

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

# OK, 1 individual that moved between plots (not counting ones that were duplicated in one or more years, which would need to resolved first) We need to fix those. Here they are:
allIndv[multPlots,1:9]

# and which ones change species
allIndv$nSp <- apply(allIndv[,Sn],1,lunique)
head(allIndv$nSp)
(multSp <- which(allIndv$nSp>1))

# what did we find?
table(allIndv$nSp)

# Oh, interesting - individuals not assigned to any species? But none assigned to more than 2! (as of 5/15/22)
probs <- which(allIndv$nSp %in% c(0,2))
length(probs)

# a bunch of them are indets seen only once in PPW1330, after the fire

# OK, 162 individuals with more than one Sp ID! We'll need to fix or exclude these. Here they are:
head(allIndv[probs,1:9])

## which individuals are missing from 2018 and present in 2013 and 2019
midNA <- function(x)
{
  if (!is.na(x[1]) & is.na(x[2]) & !is.na(x[3])) ret <- 1 else ret <- 0
  return(ret)
}
allIndv$NA18 <- apply(allIndv[,Pn[1:3]],1,midNA)
m18 <- which(allIndv$NA18==1)
allIndv$NA19 <- apply(allIndv[,Pn[2:4]],1,midNA)
m19 <- which(allIndv$NA19==1)

#make combined list
allprobs <- union(multPlots,union(multSp,union(m18,m19)))
allMults <- allIndv[allprobs,]
head(allMults)
write.csv(allMults[order(allMults$nSp,allMults$nP,decreasing = T),],'data/mult-plots-species.csv')

# write entire indvData for use in next steps
head(allIndv)
write.csv(allIndv,'data/allIndv.csv')

## Now, go back to the four data frames and match/copy these four flags, so they can be used to eliminate individuals from analysis as needed
i=1
for (i in 1:length(all.id))
{
  mt <- match(all.id[[i]]$Num,allIndv$Num)
  all.id[[i]] <- data.frame(all.id[[i]],nP=allIndv$nP[mt],nSp=allIndv$nSp[mt],NA18=allIndv$NA18[mt],NA19=allIndv$NA19[mt])
}
saveRDS(all.id,'data/allid-nodups.Rdata')

