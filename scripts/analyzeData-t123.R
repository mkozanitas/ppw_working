#analyzeData-t123 is to compare 2013, 2018 and 2019 for delayed mortality and other changes in state

# RUN prepareData and examineData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R

rm(list=ls())
library(lme4)

# Functions to convert basal area to diameter and log-diameter, and back
lba2d <- function(x) 2*sqrt((10^x)/pi)
ba2d <- function(x) 2*sqrt((x)/pi)
d2ba <- function(x) pi*(x/2)^2
d2lba <- function(x) log10(pi*(x/2)^2)

# Load fire severity data
library("RCurl")
fs <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodFireSeverity/master/data/FSextract/vegplots-54-20m-FS.csv")
head(fs)
tail(fs)
####

# load up all individual data (id) - list of 4 data.frames, one per year (133, 18, 19, 20)
all.id <- readRDS('data/allid-nodups.Rdata')
str(all.id)
length(all.id)
head(all.id[[1]])

## CONVERT Sapling diameters to adjusted values based on dbh~sadb regression
i=1
for (i in 1:length(all.id)) {
  SArows <- which(all.id[[i]]$Type=='SA')
  all.id[[i]]$BD <- NA
  all.id[[i]]$BD[SArows] <- all.id[[i]]$dbh[SArows]
  print(tail(sort(all.id[[i]]$dbh[SArows])))
  sdbh <- all.id[[i]]$dbh[SArows]
  sdbh <- sdbh * 0.64426 - 0.06439
  all.id[[i]]$dbh[SArows] <- sdbh
}

# Examine basal diameter of SAs
head(all.id[[1]])
sap13 <- all.id[[1]]
sap13 <- sap13[which(sap13$Type=='SA'),]
dim(sap13)
hist(sap13$dbh,xlim=c(0,5),breaks=c(0,1,2,3,4,5,1000))
summary(sap13$dbh)
sort(sap13$dbh[which(sap13$dbh>1)])
plot(sap13$BD,sap13$dbh,log='xy')
sap13[which(sap13$dbh>3),]
abline(0,1)
# end examine basal diameter

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

# COMNMENT OUT IF CODE ABOVE COMMENTED IN
t1 <- all.id[[1]]
t2 <- all.id[[2]]
t3 <- all.id[[3]]

names(t1)[-4] <- paste(names(t1)[-4],'.13',sep='')
names(t1)
names(t2)[-4] <- paste(names(t2)[-4],'.18',sep='')
names(t2)
names(t3)[-4] <- paste(names(t3)[-4],'.19',sep='')
names(t3)
## END TREE SELECTION SNIPPET


# And merge!
t12 <- merge(t1,t2,by = 'Num',all = T)
names(t12)

t123 <- merge(t12,t3,by = 'Num',all = T)
names(t123)
dim(t123)
head(t123)
tail(t123)

### STOPPED - deleted below this point (what was copied from t23 script and then change)
## next move from analyze data into this script using new col names

