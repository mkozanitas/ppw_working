# RUN prepareData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R
# data files are not synced by git, so prepareData does need to be run locally
rm(list=ls())
source('scripts/PW_functions_local-test.R')
source('scripts/PWfunctions_GitHub_local.R')

# read in list of 4 items, each item with the full data file for all trees in a given year
years <- c(2013,2018,2019,2020)
all.id <- readRDS('data/allidb.Rdata')

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
tail(dups)

dups
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
summary(all.id[[3]]$gCrown) # GOOD

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
# SA.patts <- c('00NANA1010','01NANA0110','10NANA0101','11NANA0101')
# TR.patts <- c('00001010','01000110','10010101','10100101','10110101','11010101','11100101','11110101')
# 
# i=2
# for (i in 2:4) {
#   all.id[[i]]$pattern <- apply(all.id[[i]][,c("Survival","bSprout","Epicormic","Apical","Dead","Live","Topkill","gCrown")],1,catVals)
#   print(table(all.id[[i]]$pattern[all.id[[i]]$Type=='SA']))
#   badSAs <- which(all.id[[i]]$Type=='SA' & !all.id[[i]]$pattern %in% SA.patts)
#   length(badSAs)
#   print(all.id[[i]][badSAs,c('Plot','Num','pattern')])
#   
#   print(table(all.id[[i]]$pattern[all.id[[i]]$Type=='TR']))
#   badTRs <- which(all.id[[i]]$Type=='TR' & !all.id[[i]]$pattern %in% TR.patts)
#   length(badTRs)
#   print(all.id[[i]][badTRs,c('Plot','Num','pattern')])
# }

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

## COMMENTED OUT - WAS A TROUBLESHOOTING STEP
# Now create some combined states, again to look for 'illegal' data combinations. Need to revisit this - some of the '2s' may suggest problems, but not sure.
# i=4
# for (i in 2:4) {
#   all.id[[i]]$Dead <- 1 - all.id[[i]]$Live
#   all.id[[i]]$DT <- all.id[[i]]$Dead + all.id[[i]]$Topkill
#   all.id[[i]]$TG <- all.id[[i]]$Topkill + all.id[[i]]$gCrown
#   all.id[[i]]$DG <- all.id[[i]]$Dead + all.id[[i]]$gCrown
#   all.id[[i]]$TB <- 0
#   all.id[[i]]$TB[which(all.id[[i]]$bSprout==1 & all.id[[i]]$Topkill==1)] <- 1
#   print(c(years[i],'Dead'))
#   print(table(all.id[[i]]$Dead))
#   print(c(''))
#   print(c(years[i],'Topkill'))
#   print(table(all.id[[i]]$Topkill))
#   print(c(''))
#   print(c(years[i],'gCrown'))
#   print(table(all.id[[i]]$gCrown))
#   print(c(years[i],'DT'))
#   print(table(all.id[[i]]$DT))
#   print(c(years[i],'TG'))
#   print(table(all.id[[i]]$TG))
#   print(c(years[i],'DG'))
#   print(table(all.id[[i]]$DG))
# }

## CREATE FOUR FATES
i=4
for (i in 2:4) {
  all.id[[i]]$DN <- 0
  all.id[[i]]$DR <- 0
  all.id[[i]]$LN <- 0
  all.id[[i]]$LR <- 0
  
  # all.id[[i]]$DN[which(all.id[[i]]$Topkill==0 & all.id[[i]]$bSprout==1)] <- 0
  # all.id[[i]]$DR[which(all.id[[i]]$Topkill==0 & all.id[[i]]$bSprout==0)] <- 0
  # all.id[[i]]$LN[which(all.id[[i]]$Topkill==1 & all.id[[i]]$bSprout==1)] <- 0
  # all.id[[i]]$LR[which(all.id[[i]]$Topkill==1 & all.id[[i]]$bSprout==0)] <- 0
  
  all.id[[i]]$DN[which(all.id[[i]]$Topkill==1 & all.id[[i]]$bSprout==0)] <- 1
  all.id[[i]]$DR[which(all.id[[i]]$Topkill==1 & all.id[[i]]$bSprout==1)] <- 1
  all.id[[i]]$LN[which(all.id[[i]]$Topkill==0 & all.id[[i]]$bSprout==0)] <- 1
  all.id[[i]]$LR[which(all.id[[i]]$Topkill==0 & all.id[[i]]$bSprout==1)] <- 1
}
i=2
nrow(all.id[[i]])
table(all.id[[i]]$DN,useNA = 'always')
table(all.id[[i]]$DR,useNA = 'always')
table(all.id[[i]]$LN,useNA = 'always')
table(all.id[[i]]$LR,useNA = 'always')

i=3
for (i in 2:4) {
  print(table(all.id[[i]][,c('Topkill','gCrown')]))
}


## SKIP AS ALL DATA CLEAN AS OF 5/23/23
if (FALSE) {
  for (i in 2:4) {
    print(table(all.id[[i]][,c('DN','DR')]))
    print(table(all.id[[i]][,c('DN','LN')]))
    print(table(all.id[[i]][,c('DN','LR')]))
    print(table(all.id[[i]][,c('DR','LN')]))
    print(table(all.id[[i]][,c('DR','LR')]))
    print(table(all.id[[i]][,c('LN','LR')]))
  }
  
  # Every tree should be one of these fates
  # SUCCESS!! as of 5/31/23
  for (i in 2:4) {
    all.id[[i]]$nFates <- apply(all.id[[i]][,c('DN','DR','LN','LR')],1,sum)
    print(table(all.id[[i]]$nFates,useNA = 'always'))
  }
  
  tfail <- which(all.id[[i]]$nFates==0)
  length(tfail)
  head(all.id[[i]][tfail,])
}

## Create fates variable with DN, DR, LN, LR as states
for (i in 2:4) {
  all.id[[i]]$fate <- NA
  all.id[[i]]$fate[which(all.id[[i]]$DN==1)] <- 'DN'
  all.id[[i]]$fate[which(all.id[[i]]$DR==1)] <- 'DR'
  all.id[[i]]$fate[which(all.id[[i]]$LN==1)] <- 'LN'
  all.id[[i]]$fate[which(all.id[[i]]$LR==1)] <- 'LR'
  print(table(all.id[[i]]$fate,useNA='always'))
}

# check that our previous combined fates align with these new ones
for (i in 2:4) print(table(all.id[[i]]$Live,all.id[[i]]$fate,useNA='always'))

for (i in 2:4) print(table(all.id[[i]]$gCrown,all.id[[i]]$fate,useNA='always'))

for (i in 2:4) print(table(all.id[[i]]$Topkill,all.id[[i]]$fate,useNA='always'))

for (i in 2:4) {
  all.id[[i]]$Resprout <- NA
  all.id[[i]]$Resprout[which(all.id[[i]]$fate %in% c('DR','LR'))] <- 1
  all.id[[i]]$Resprout[which(all.id[[i]]$fate %in% c('DN','LN'))] <- 0
  print(table(all.id[[i]]$Resprout,all.id[[i]]$fate,useNA='always'))
}

# a few more problem plants!
i=2
all.id[[i]]$Num[which(all.id[[i]]$fate %in% c('LN','LR') & all.id[[i]]$gCrown==0)]

i=3
all.id[[i]]$Num[which(all.id[[i]]$gCrown==-Inf)]
all.id[[i]]$Num[which(all.id[[i]]$fate %in% c('LN') & all.id[[i]]$gCrown==0)]


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

# load up all individual data (id) - list of 4 data.frames, one per year (133, 18, 19, 20)
str(all.id)
length(all.id)
head(all.id[[1]])
head(all.id[[1]][all.id[[1]]$Type=='SA',])

# how many trees have BD but not DBH
nodbh <- which(!is.na(all.id[[2]]$SA.BD_cm) & is.na(all.id[[2]]$DBH_cm))
table(all.id[[2]]$Type[nodbh])
length(nodbh)

## CREATE A CALCULATED BASAL DIAMETER AT 10 CM FOR TREES (b10)
if (TRUE) {
  i=1
  for (i in 1:length(all.id)) {
    TRrows <- which(all.id[[i]]$Type=='TR')

    all.id[[i]]$d10 <- all.id[[i]]$SA.BD_cm
    # summary(all.id[[i]]$d10,useNA='always')
    
    # from DBH-SADB.R script
    # D10 = DBH.cm * 1.176 + 1.070
    all.id[[i]]$d10[TRrows] <- all.id[[i]]$dbh[TRrows] * 1.176 + 1.07
    summary(all.id[[i]]$d10,useNA='always')
}
  
  # Examine basal diameter of SAs
  head(all.id[[1]])
  sap13 <- all.id[[1]]
  sap13 <- sap13[which(sap13$Type=='SA'),]
  dim(sap13)
  hist(sap13$dbh)
  summary(sap13$dbh)
  length(which(sap13$dbh<0.01))
  nrow(sap13)
  hist(sap13$SA.Height_cm)
  plot(sap13$dbh,sap13$SA.Height_cm,xlim=c(-1,2))
  
  sort(sap13$dbh[which(sap13$dbh>1)])
  plot(sap13$SA.BD_cm,sap13$dbh,log='')
  sap13[which(sap13$dbh>3),]
  abline(0,1)
  # end examine basal diameter
}

spNames <- read.csv('data/all-spp-names.csv')
head(spNames)
names(spNames)[which(names(spNames)=='x')] <- 'spName'
#allIndv <- readRDS('data/allIndv.Rdata')
allIndv <- read.csv('data/allIndv.csv')
head(allIndv)

#### First round of analysis of post-fire states, 2013-2018
## get percent survival by species and type

# THIS SNIPPET WOULD IDENTIFY TREES THAT CHANGE PLOT OR ID AND SHOULD BE ELIMINATED - FOR NOW ANALYZING EVERYTHING
# take any individuals where plot and species match for the first two years
nrow(allIndv)
#plot.ok <- NA
#spec.ok <- NA
#ps.ok <- intersect(plot.ok,spec.ok)
#length(ps.ok)
#head(ps.ok)
#tail(ps.ok)

# include trees from new plots
# ps.ok <- c(ps.ok,allIndv$Num[which(allIndv$P18 %in% c('PPW1851','PPW1852','PPW1853','PPW1854'))])
# length(ps.ok)
# head(ps.ok)
# tail(ps.ok)

# use this if we want to subset to valid trees - inactivated for now
# t1 <- all.id[[1]][which(all.id[[1]]$Num %in% ps.ok),]
# nrow(t1)
# t2 <- all.id[[2]][which(all.id[[2]]$Num %in% ps.ok),]
# nrow(t2)
# t3 <- all.id[[3]][which(all.id[[3]]$Num %in% ps.ok),]
# nrow(t3)
## END TREE SELECTION SNIPPET

# MERGE YEARS!!
t1 <- all.id[[1]]
t2 <- all.id[[2]]
t3 <- all.id[[3]]
t4 <- all.id[[4]]

names(t1)
names(t1)[-4] <- paste(names(t1)[-4],'.13',sep='')
names(t1)
names(t2)[-4] <- paste(names(t2)[-4],'.18',sep='')
names(t2)
names(t3)[-4] <- paste(names(t3)[-4],'.19',sep='')
names(t3)
names(t4)[-4] <- paste(names(t4)[-4],'.20',sep='')
names(t4)

# And merge!
t12 <- merge(t1,t2,by = 'Num',all = T)
names(t12)

t123 <- merge(t12,t3,by = 'Num',all = T)
names(t123)
dim(t123)
head(t123)
tail(t123)

tAll <- merge(t123,t4,by = 'Num',all = T)
names(tAll)
dim(tAll)
head(tAll)
tail(tAll)

rm('t12')
rm('t123')

# All four years are now merged!!!!! ###

# create 'proxy 2013 data'
# rownums for new 2018 individuals from new plots
table(tAll$Plot.13)
table(tAll$Plot.18)

# three subsets of new individuals
# we assume that all newIndvs were present just before the fire, as we either tagged them alive and recovering or dead; all of these should be included in estimates of fates
# newIndvs: new recruits, and recruited and dead, and newplots
newIndvs <- which(is.na(tAll$Plot.13))
length(newIndvs)
table(tAll$Type.18[newIndvs])
summary(tAll$Num[newIndvs])
length(which(tAll$Num[newIndvs]>99000))

# newly recruited and dead; subset of newIndvs
n99 <- which(tAll$Num>99000)
length(n99)
table(tAll$Type.18[n99])

# new plots; subset of newIndvs
newPlot <- which(tAll$Plot.18 %in% c('PPW1851','PPW1852','PPW1853','PPW1854'))
length(newPlot)
table(tAll$Type.18[newPlot])

# New Trees - either not measured in 2013 or new recruits (or other error?)
newTrees <- which(!tAll$Num %in% tAll$Num[n99] & !tAll$Num %in% tAll$Num[newPlot] & tAll$Type.18=='TR' & is.na(tAll$Plot.13) & tAll$Num %% 1==0)
length(newTrees)
sort(tAll$Num[newTrees])
table(tAll$Plot.18[newTrees])

#tAll[newTrees[15],]
# looks like 12 trees with Num < 5500 that were tagged and no data collected in 2013 - then there are a few new trees that may have recruited from <sapling to tree stage

length(intersect(newIndvs,n99))
length(intersect(newIndvs,newPlot))
length(intersect(n99,newPlot))

summary(tAll$dbh.18[newPlot])
summary(tAll$dbh.18[n99])

# Now start filling in pre-fire data from 2018 data, for completeness, for all new individuals
#create13Data <- c(newPlot,newTrees[1:12])
tAll$Plot.13[newIndvs] <- tAll$Plot.18[newIndvs]
tAll$Quad.13[newIndvs] <- tAll$Quad.18[newIndvs]
tAll$Type.13[newIndvs] <- tAll$Type.18[newIndvs]
tAll$Species.13[newIndvs] <- tAll$Species.18[newIndvs]

tAll$Dead.13[newIndvs] <- 0
tAll$Live.13[newIndvs] <- 1
tAll$gCrown.13[newIndvs] <- 1
tAll$dbh.13[newIndvs] <- tAll$dbh.18[newIndvs]

# ignore these individuals for basal area growth - commented out because we aren't creating proxy 2013 diameter data now
tAll$UseForBAGrowth <- T
tAll$UseForBAGrowth[newIndvs] <- F

tAll$TreeNum <- floor(tAll$Num)
tAll$fPlot <- as.factor(tAll$Plot.18)

# Load plot characteristics
northness <- function(asp,slp,units='deg') {
  if (units=='deg') {
    asp <- 2*pi*asp/360
    slp <- 2*pi*slp/360
  }
  cos(asp)*sin(slp)
}

eastness <- function(asp,slp,units='deg') {
  if (units=='deg') {
    asp <- 2*pi*asp/360
    slp <- 2*pi*slp/360
  }
  sin(asp)*sin(slp)
}

plotInfo <- get.plot.info()
plotInfo
plotInfo$northness <- northness(plotInfo$Aspect,plotInfo$Slope)
plotInfo$eastness <- eastness(plotInfo$Aspect,plotInfo$Slope)
###

# transfer northness to tAll
p2t <- match(tAll$Plot.18,plotInfo$Plot)
tAll$northness <- plotInfo$northness[p2t]
tAll$eastness <- plotInfo$eastness[p2t]

write.csv(tAll,'data/tAll.csv')

