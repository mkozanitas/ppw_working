# RUN prepareData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R
# data files are not synced by git, so prepareData does need to be run locally
rm(list=ls())
source('scripts/11.PW_functions_local-test.R')
source('scripts/12.PW_functions_GitHub_local.R')

# read in list of 4 items, each item with the full data file for all trees in a given year. We are using individual points for analysis in this paper (hence 'allidb.Rdata')
years <- c(2013,2018,2019,2020)
all.id <- readRDS('data/allidb.Rdata')

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
  #if (length(dN)>0) all.id[[i]] <- all.id[[i]][-which(all.id[[i]]$Num %in% dN),]
}

dups[[1]]
dups <- do.call(rbind,dups)
dim(dups)
head(dups)
tail(dups)

dups
write.csv(dups,'data/duplicates.csv')

## Check on 12/24/24 showed two dups in 2013. For the first, the first line matches the 2018 data (Species) and for the second, they seem identical except the first line has size data and the 2nd doesn't. So drop the second in each case
all.id[[1]] <- all.id[[1]][-c(2973,2974),]

# check that the maximum numbers from each year of survey don't have typos or bad values
for (i in 1:4) print(tail(sort(all.id[[i]]$Num)))

# Pull out all individuals with NA for bSprout - this should always be filled out. As of 5/15/22, no problems identified at this step!!
# 9/24/24 - Melina: There are 944 NA for bSprout in 2013. Is that okay?
bNA <- all.id[[1]][which(is.na(all.id[[1]]$bSprout)),]
for (i in 2:4) bNA <- rbind(bNA,all.id[[i]][which(is.na(all.id[[i]]$bSprout)),])
dim(bNA)
head(bNA)
table(bNA$Year)
write.csv(bNA,'data/bSprout-NAs.csv')

## Seems there are problems with gCrown - checking here
table(all.id[[1]]$gCrown,useNA='always') # GOOD
table(all.id[[2]]$gCrown,useNA='always') # a bunch of NAs
table(all.id[[3]]$gCrown,useNA='always') # a bunch of NAs

summary(all.id[[4]]$gCrown) # Character! Let's make it numeric
all.id[[4]]$gCrown <- as.numeric(all.id[[4]]$gCrown)
table(all.id[[4]]$gCrown,useNA='always') # a bunch of NAs

## Assign Dead (Live=0), and fix gCrown NAs
all.id[[1]]$Dead <- 1 - all.id[[1]]$Live

i=4
for (i in 2:4)
{
  all.id[[i]]$Dead <- 1 - all.id[[i]]$Live
  all.id[[i]]$Epi.Api <- 0
  all.id[[i]]$Epi.Api[which(all.id[[i]]$Epicormic==1 | all.id[[i]]$Apical==1)] <- 1
  
  # fix gCrown NAs
  if (FALSE) {
    head(all.id[[i]][which(all.id[[i]]$Live==1 & is.na(all.id[[i]]$gCrown)),c('Plot','Num','Survival','bSprout','Topkill','Live','gCrown')])
    table(all.id[[i]]$Survival,all.id[[i]]$Topkill,useNA='always')
  }
  
  all.id[[i]]$gCrown[which(all.id[[i]]$Live==1)] <- 
    all.id[[i]]$Survival[which(all.id[[i]]$Live==1)]
  all.id[[i]]$gCrown[which(all.id[[i]]$Live==0)] <- 0
  table(all.id[[i]]$gCrown,useNA='always')
} 


## TROUBLESHOOTING COMBINATIONS, TO COME UP WITH GCROWN FIX THAT'S IN LOOP ABOVE
if (FALSE) {
  i=2 
  table(all.id[[i]]$Epicormic,all.id[[i]]$Apical,all.id[[i]]$Epi.Api,useNA='always')
  table(all.id[[i]]$Survival,all.id[[i]]$Epi.Api,useNA='always')
  table(all.id[[i]]$gCrown,all.id[[i]]$Epi.Api,useNA='always')
  
  print(table(all.id[[i]]$Live,all.id[[i]]$Survival,useNA='always'))
  print(table(all.id[[i]]$Survival,all.id[[i]]$gCrown,useNA='always'))
  print(table(all.id[[i]]$gCrown,all.id[[i]]$Topkill,all.id[[i]]$Live,useNA='always'))
  print(table(all.id[[i]]$gCrown,all.id[[i]]$Survival,all.id[[i]]$bSprout,useNA='always'))
  
  
  # In post-fire years, check individuals scored as any combination of DEAD & TOPKILL, DEAD & GREEN, TOPKILL & GREEN
  i=2
  names(all.id[[i]])
  table(all.id[[i]]$Dead,all.id[[i]]$Live,useNA='always')
  table(all.id[[i]]$gCrown,all.id[[i]]$Live,useNA='always')
  
  # WHO ARE THESE INDIVIDUALS WITH NA FOR gCrown and live or dead
  head(all.id[[i]][which(is.na(all.id[[i]]$gCrown) & all.id[[i]]$Live==0),])
  
  table(all.id[[i]]$Topkill,all.id[[i]]$Live,useNA='always')
  table(all.id[[i]]$bSprout,all.id[[i]]$Live,useNA='always')
  
  # catenate values to see patterns
  catVals <- function(x) {
    res <- c()
    for (i in 1:length(x)) res <- paste(res,x[i],sep='')
    return(res)
  }
  
  ## ONLY DO THIS FOR 2018 AND BEYOND (i in 2:4)
  # These are patterns of values for 8 fields in the data (see below) which represent the 'legal' combinations. Any tree that doesn't follow this pattern suggests either a data ehtry or a coding error re
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
}

## CREATE FOUR FATES
i=4
for (i in 2:4) {
  all.id[[i]]$DN <- 0
  all.id[[i]]$DN[which(all.id[[i]]$Topkill==1 & all.id[[i]]$bSprout==0)] <- 1
  
  all.id[[i]]$DR <- 0
  all.id[[i]]$DR[which(all.id[[i]]$Topkill==1 & all.id[[i]]$bSprout==1)] <- 1
  
  all.id[[i]]$LN <- 0
  all.id[[i]]$LN[which(all.id[[i]]$Topkill==0 & all.id[[i]]$bSprout==0)] <- 1
  
  all.id[[i]]$LR <- 0
  all.id[[i]]$LR[which(all.id[[i]]$Topkill==0 & all.id[[i]]$bSprout==1)] <- 1
  
  all.id[[i]]$fate4 <- NA
  all.id[[i]]$fate4[which(all.id[[i]]$DN==1)] <- 'DN'
  all.id[[i]]$fate4[which(all.id[[i]]$DR==1)] <- 'DR'
  all.id[[i]]$fate4[which(all.id[[i]]$LN==1)] <- 'LN'
  all.id[[i]]$fate4[which(all.id[[i]]$LR==1)] <- 'LR'
  
  all.id[[i]]$fate3 <- all.id[[i]]$fate4
  all.id[[i]]$fate3[which(all.id[[i]]$fate4 %in% c('LN','LR'))] <- 'GC'
  
  all.id[[i]]$Resprout <- NA
  all.id[[i]]$Resprout[which(all.id[[i]]$fate4 %in% c('DR','LR'))] <- 1
  all.id[[i]]$Resprout[which(all.id[[i]]$fate4 %in% c('DN','LN'))] <- 0
}

## SKIP AS ALL DATA CLEAN AS OF 5/23/23
if (FALSE) {
  i=2
  for (i in 2:4) 
  {
    nrow(all.id[[i]])
    table(all.id[[i]]$DN,useNA = 'always')
    table(all.id[[i]]$DR,useNA = 'always')
    table(all.id[[i]]$LN,useNA = 'always')
    table(all.id[[i]]$LR,useNA = 'always')
    table(all.id[[i]]$fate4,useNA = 'always')
    table(all.id[[i]]$fate3,useNA = 'always')
    
    print(table(all.id[[i]][,c('DN','DR')]))
    print(table(all.id[[i]][,c('DN','LN')]))
    print(table(all.id[[i]][,c('DN','LR')]))
    print(table(all.id[[i]][,c('DR','LN')]))
    print(table(all.id[[i]][,c('DR','LR')]))
    print(table(all.id[[i]][,c('LN','LR')]))
    
    print(table(all.id[[i]]$Live,all.id[[i]]$fate4,useNA='always'))
    print(table(all.id[[i]]$gCrown,all.id[[i]]$fate4,useNA='always'))
    print(table(all.id[[i]]$Topkill,all.id[[i]]$fate4,useNA='always'))
    print(table(all.id[[i]]$Resprout,all.id[[i]]$fate4,useNA='always'))
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

# 12/24/24 - just one problem lefgt
allIndv[3726,]

if (FALSE) ### All of these have been fixed
{
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
}

# write entire indvData for use in next steps
head(allIndv)
write.csv(allIndv,'data/allIndv.csv')

## Now, go back to the four data frames and match/copy these four flags, so they can be used to eliminate individuals from analysis as needed
i=1
for (i in 1:length(all.id))
{
  mt <- match(all.id[[i]]$Num,allIndv$Num)
  all.id[[i]] <- data.frame(all.id[[i]],nP=allIndv$nP[mt],nSp=allIndv$nSp[mt])
}
saveRDS(all.id,'data/allid-nodups.rds')

# how many trees have BD but not DBH
nodbh <- which(!is.na(all.id[[1]]$SA.BD_cm) & is.na(all.id[[1]]$DBH_cm))
table(all.id[[1]]$Type[nodbh])
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
    all.id[[i]]$d10[TRrows] <- all.id[[i]]$DBH_cm[TRrows] * 1.176 + 1.07
    summary(all.id[[i]]$d10,useNA='always')
  }
  
  # Examine basal diameter of SAs
  if (FALSE) {
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
    
    tr13 <- all.id[[1]]
    tr13 <- tr13[which(tr13$Type=='TR'),]
    names(tr13)
    plot(tr13$DBH_cm,tr13$dbh,log='')
    # end examine basal diameter
  }
}

spNames <- read.csv('data/all-spp-names.csv')
head(spNames)
tail(spNames)
names(spNames)[which(names(spNames)=='x')] <- 'spName'
#allIndv <- readRDS('data/allIndv.Rdata')
allIndv <- read.csv('data/allIndv.csv')
head(allIndv)

# MERGE YEARS!!
t0 <- all.id[[1]]
t1 <- all.id[[1]]
t2 <- all.id[[2]]
t3 <- all.id[[3]]
t4 <- all.id[[4]]

# how many individuals in each year, before merging
nrow(t0)
nrow(t1)
nrow(t2)
nrow(t3)
nrow(t4)

names(t0)
names(t0)[-4] <- paste(names(t1)[-4],'.13',sep='')
names(t1)
names(t1)[-4] <- paste(names(t1)[-4],'.17',sep='')
names(t1)
names(t2)[-4] <- paste(names(t2)[-4],'.18',sep='')
names(t2)
names(t3)[-4] <- paste(names(t3)[-4],'.19',sep='')
names(t3)
names(t4)[-4] <- paste(names(t4)[-4],'.20',sep='')
names(t4)

# And merge!
t01 <- merge(t0,t1,by = 'Num',all = T)
names(t01)

t12 <- merge(t01,t2,by = 'Num',all = T)
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

rm('t01')
rm('t12')
rm('t123')

# create master Plot variable
tAll$Plot <- tAll$Plot.13
tAll$Plot[which(is.na(tAll$Plot))] <- tAll$Plot.18[which(is.na(tAll$Plot))]
tAll$Plot[which(is.na(tAll$Plot))] <- tAll$Plot.19[which(is.na(tAll$Plot))]
tAll$Plot[which(is.na(tAll$Plot))] <- tAll$Plot.20[which(is.na(tAll$Plot))]
substr(tAll$Plot,4,5) <- '13'
table(tAll$Plot,useNA='always')

# create master Species variable - note there was one individual with a species problem, and this assigns it to first value
tAll$Species <- tAll$Species.13
tAll$Species[which(is.na(tAll$Species))] <- tAll$Species.18[which(is.na(tAll$Species))]
tAll$Species[which(is.na(tAll$Species))] <- tAll$Species.19[which(is.na(tAll$Species))]
tAll$Species[which(is.na(tAll$Species))] <- tAll$Species.20[which(is.na(tAll$Species))]
table(tAll$Species,useNA='always')

# Create TreeNum variable identify individuals, removing point numbers
tAll$TreeNum <- floor(tAll$Num)
head(table(tAll$TreeNum))

# create 'proxy 2017 data'

# identify several subsets of individuals
# present and measured in 2013
# inferred pre-fire in old plots, or tagged and not recorded in 2013, and dead (99000s and some lower numbers)
# inferred pre-fire, old plots, and alive (new tags)
# tagged in 2013 and not written down in data
# new plot, alive after the fire

tAll$Cat17 <- NA

# tagged and recorded in data sheets in 2013
tAll$Cat17[which(!is.na(tAll$Plot.13))] <- 'Tag13'

# not recorded in 2013 (some tagged!) and dead, mostly 99000s
tAll$Cat17[which(is.na(tAll$Live.13) & tAll$Live.18==0)] <- '99s.old'

# not present in 2013, found alive in 18; or tagged but not recorded in 2013, found alive in 18
tAll$Cat17[which(is.na(tAll$Live.13) & tAll$Live.18==1)] <- 'new17'

# new plots
tAll$Cat17[which(tAll$Plot.18 %in% c('PPW1851','PPW1852','PPW1853','PPW1854'))] <- 'NewPlot'

tAll$Cat17[which(is.na(tAll$Live.13) & is.na(tAll$Live.18) & tAll$Live.19==1)] <- 'new19'

tAll$Cat17[which(is.na(tAll$Live.13) & is.na(tAll$Live.18) & is.na(tAll$Live.19) & tAll$Live.20==1)] <- 'new20'

# summarize 2017 category data
table(tAll$Cat17,useNA='always')

# Not in any category - Check one at a time - DONE
#tAll[which(is.na(tAll$Cat17))[8],]

# create flag to exclude problematic stems
tAll$ExcStem <- 0
tAll$ExcStem[which(is.na(tAll$Cat17))[c(2:5,8)]] <- 1

# one QUEAGR without a num
tAll[tAll$Plot=='PPW1339',c('Species.13','Num')]

# what plot were 99s found in
table(tAll$Plot.18[which(tAll$Num>10000)])
table(tAll$Plot.18[which(tAll$Cat17 == '99s.old')])
table(tAll$Plot.18[which(tAll$Cat17 == 'NewPlot')])

# three subsets of new individuals
# we assume that all newIndvs were present just before the fire, as we either tagged them alive and recovering or dead; all of these should be included in estimates of fates
# newIndvs: new recruits, and recruited and dead, and newplots
newIndvs <- which(is.na(tAll$Plot.13) & !is.na(tAll$Plot.18))

# newIndvs matches the three categories above
table(tAll$Cat17[newIndvs])
table(tAll$Cat17,useNA='always')
all(sort(newIndvs),sort(which(tAll$Cat17 %in% c('NewPlot','new17','99s.old'))))

length(newIndvs) # matches Cat17 which weren't tagged in 2013

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

# flag to identify plants that really are post-fire recruits, and not overlooked plants that were there before the fire
tAll$PFRecruit <- 0

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

# how many new plants added each year
table(tAll$Type.13,tAll$Type.18,useNA = 'always')
table(tAll$Type.18,tAll$Type.19,useNA = 'always')
table(tAll$Type.19,tAll$Type.20,useNA = 'always')

# where are the new plants
new18 <- which(is.na(tAll$Type.13) & !is.na(tAll$Type.18))
length(new18)
new19 <- which(is.na(tAll$Type.18) & !is.na(tAll$Type.19))
length(new19)
new20 <- which(is.na(tAll$Type.19) & !is.na(tAll$Type.20))
length(new20)

(p18 <- table(tAll$Plot[new18],useNA='always'))
sum(p18)

# how many new plants in previously surveyed plots - found in 44 of 50 plots:
sum(p18[1:44])

table(tAll$Plot[new19],useNA='always')
table(tAll$Plot[new20],useNA='always')

table(tAll$Plot[new18],tAll$Species[new18])
table(tAll$Plot[new19],tAll$Species[new19])
table(tAll$Plot[new20],tAll$Species[new20])

hist(tAll$DBH_cm.18[new18])
hist(tAll$DBH_cm.19[new19])
hist(tAll$DBH_cm.20[new20])

# output suspicious transitions- can switch Type and year here to look at other combos
# MELINA_CHECK
s1 <- which(tAll$Type.13=='TS' & tAll$Type.18=='TR')
tAll[s1,c('Plot.13','Num','Species.13')]

s1 <- which(tAll$Type.18=='TS' & tAll$Type.19=='TR')
tAll[s1,c('Plot.13','Num')]

s1 <- which(tAll$Type.19=='TS' & tAll$Type.20=='TR')
tAll[s1,c('Plot.13','Num')]

new19sap <- which(is.na(tAll$Type.18) & !is.na(tAll$Type.19))
tAll[new19sap,c('Plot.19','Num')]
table(tAll$Plot.19[new19sap])
table(tAll$Type.19[new19sap])


# Now start filling in pre-fire data from 2018 data, to complete 2017 proxy data, for all new individuals
#create13Data <- c(newPlot,newTrees[1:12])
table(tAll$Plot.17)
tAll$Plot.17[newIndvs] <- tAll$Plot.18[newIndvs]
tAll$Quad.17[newIndvs] <- tAll$Quad.18[newIndvs]
tAll$Type.17[newIndvs] <- tAll$Type.18[newIndvs]

# changing type to TR if the new size shows it recruited up
tc <- which(tAll$Type.18=='TR' & tAll$Type.17=='SA')
tAll$Type.17[tc] <- tAll$Type.18[tc]

tAll$Species.17[newIndvs] <- tAll$Species.18[newIndvs]

# all assumed to be alive right before fire. Introduces small error due to mortality between 2013 and 2017
tAll$Dead.17[newIndvs] <- 0
tAll$Live.17[newIndvs] <- 1
tAll$gCrown.17[newIndvs] <- 1

# Now, what to do about size! Create unified d10 variable in 2017 for models. And transfer SA and tree DBH as well
length(which(!is.na(tAll$d10.13)))
length(which(!is.na(tAll$Plot.13)))

# start by assigning 2018 dbh to 2017 proxy data
tAll$d10.17 <- tAll$d10.18
tAll$DBH_cm.17 <- tAll$DBH_cm.18
tAll$SA.BD_cm.17 <- tAll$SA.BD_cm.18

# then use 2013 if 2018 was missing, or if diameter shrunk - 995 trees
rsel <- which(is.na(tAll$d10.17) & !is.na(tAll$d10.13))
shr <- which(tAll$d10.18<tAll$d10.13)
rsel <- union(rsel,shr)
tAll$d10.17[rsel] <- tAll$d10.13[rsel]

rsel <- which(is.na(tAll$DBH_cm.17) & !is.na(tAll$DBH_cm.13))
shr <- which(tAll$DBH_cm.18<tAll$DBH_cm.13)
rsel <- union(rsel,shr)
tAll$DBH_cm.17[rsel] <- tAll$DBH_cm.13[rsel]

rsel <- which(is.na(tAll$SA.BD_cm.17) & !is.na(tAll$SA.BD_cm.13))
shr <- which(tAll$SA.BD_cm.18 < tAll$SA.BD_cm.13)
rsel <- union(rsel,shr)
tAll$SA.BD_cm.17[rsel] <- tAll$SA.BD_cm.13[rsel]

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

## Write results to master data file!!
tAll$Survey <- 'Plot'
table(tAll$Plot)
names(tAll)
tAll$Plot.Orig <- tAll$Plot
substr(tAll$Plot,4,5) <- '13'
write.csv(tAll,'data/tAllPlot.csv')

# now add hectares for expanded data set
htAll <- readRDS('data/hectares.rds')
names(htAll)
dim(htAll)

inBoth <- intersect(tAll$Num,htAll$Num)
inBoth

# remove from htAll
# htAll <- htAll[-which(htAll$Num %in% inBoth),]
dim(htAll)

names(tAll)
names(htAll)
cbind(names(htAll),names(htAll) %in% names(tAll)
)

dim(tAll)
dim(htAll)
tAllh <- merge(tAll,htAll,by=names(htAll)[which(names(htAll) %in% names(tAll))],all=T)
dim(tAllh)
tAllh <- tAllh[-which(tAllh$ExcStem==1),]
dim(tAllh)

names(tAllh)
tail(tAllh)
table(tAllh$ExcStem,useNA='always')
table(tAllh$Plot,useNA='always')
table(tAllh$Plot.Orig.18,useNA='always')
table(tAllh$Live.17,useNA='always')
table(tAllh$Live.17,tAllh$Survey,useNA='always')

write.csv(tAllh,'data/tAllh.csv')
