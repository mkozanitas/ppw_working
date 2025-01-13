# Extract plot data for 2013, 2018, 2019, 2020 and prepare for analysis
rm(list=ls())

source('scripts/11.PW_functions_local-test.R')
source('scripts/13.hectare_scripts.r')

## FOR NOW CODE TO EXPORT DATA WITH BRANCHES ALL COMMENTED OUT

# Load the per-year data (without aggregating branches), then set to branches=F to collapse points
id13 <- get.indv.data(year = 2013,branches=F,stump=F,orig.dead=F)
dim(id13)
head(id13)
names(id13)
head(sort(id13$Num))

#to see all individual branches...points not collapsed use branches=T
id13b <- get.indv.data(year = 2013,branches=T,stump=F,orig.dead=F)
dim(id13b)
head(id13b)
head(sort(id13b$Num))
id13b[id13b$Num==1671,]

# sloppy coding to be able to turn a code snippet on and off
if (FALSE) 
{
  id13f <- get.indv.data(year = 2013, stump=T, orig.dead=T, survival=T, bsprout=T, epicormic=T, apical=T, canopy=T, bsprout.height=T, bsprout.count=T, tag.pulled=T, keep.999=T, branches=F)
  dim(id13f)
  head(id13f)
}

#indv.data.2013.B <- get.indv.data(year = 2013,branches=T)
#dim(indv.data.2013.B)
#head(indv.data.2013.B)

id18 <- get.indv.data(year = 2018,branches=F,keep.999 = T)
dim(id18)
head(id18)

id18b <- get.indv.data(year = 2018,branches=T,keep.999 = T)
dim(id18b)
head(id18b)

## run 50.troubleshoot interactively
####
id19 <- get.indv.data(year = 2019,branches=F)
dim(id19)
head(id19)

id19b <- get.indv.data(year = 2019,branches=T)
dim(id19b)
head(id19b)

plot.list.2020 <- get.plot(2020)
# now remove plots not sampled
plot.list.2020 <- plot.list.2020[-c(5,6,20,25,28,30,40:42,46,49,51:54)]
plot.list.2020

id20 <- get.indv.data(year = 2020,branches=F,plot.list=plot.list.2020)
dim(id20)
head(id20)

id20b <- get.indv.data(year = 2020,branches=T,plot.list=plot.list.2020)
dim(id20b)
head(id20b)

## write out all the data to an RData file, so it can be read and analyzed without going back to github
all.id <- list(id13,id18,id19,id20)
all.idb <- list(id13b,id18b,id19b,id20b)
saveRDS(all.id,'data/allid.rds')
saveRDS(all.idb,'data/allidb.rds')

## get hectare tree data
ht16 <- read.csv('data/hectares/hectare_trees_all_2016.csv')
dim(ht16)
head(ht16)
head(ht16$Num)
table(ht16$Plot)
plot.trees <- grep('from offset',ht16$Accuracy)
length(plot.trees)
ht16 <- ht16[-plot.trees,]

ht18 <- read.csv('data/hectares/Hectare_trees_all_2018.csv')
dim(ht18)
head(ht18)
names(ht18)
head(ht18$Num)
table(ht18$Plot)

names(ht16)[-4] <- paste(names(ht16)[-4],'.15',sep='')
names(ht18)[-4] <- paste(names(ht18)[-4],'.18',sep='')

dim(ht16)
dim(ht18)
htAll <- merge(ht16,ht18,by='Num',all = T)
dim(htAll)
names(htAll)

misPlot <- (which(htAll$Plot.15!=htAll$Plot.18))
misPlot
htAll[misPlot,c('Num','Plot.15','Subplot.15','Species.15','Plot.18','Subplot.18','Species.18')]
# remove these
htAll <- htAll[-misPlot,]

misSpp <- which(htAll$Species.15 != htAll$Species.18)
misSpp
table(htAll$Species.15[misSpp],htAll$Species.18[misSpp],useNA='always')

# how many of the misSpp are due to senescence
table(htAll$Senesced.18,useNA='always')
notSen <- which(is.na(htAll$Senesced.18))
length(notSen)

misSp2 <- which(htAll$Species.15[notSen] != htAll$Species.18[notSen])
misSp2
table(htAll$Species.15[misSp2],htAll$Species.18[misSp2],useNA='always')
# All of them!

# so, for analysis, just rely on Species.15

# how well do sizes match up?
plot(htAll$DBH_cm.15,htAll$DBH_cm.18,log='xy')
htAll[which.max(htAll$DBH_cm.18),]

# definitely looks like some errors, but I'll just rely on DBH.15, unless missing (mostly 2021 new superplot, presumably)
htAll$DBH_cm.17 <- htAll$DBH_cm.15
length(which(is.na(htAll$DBH_cm.17)))
table(htAll$Plot.18[which(is.na(htAll$DBH_cm.17))])

htAll$DBH_cm.17[which(is.na(htAll$DBH_cm.17))] <- htAll$DBH_cm.18[which(is.na(htAll$DBH_cm.17))]
length(which(is.na(htAll$DBH_cm.17)))

htAll$Species <- htAll$Species.15
length(which(is.na(htAll$Species)))
htAll$Species[which(is.na(htAll$Species))] <- htAll$Species.18[which(is.na(htAll$Species))]

htAll$Plot <- htAll$Plot.15
length(which(is.na(htAll$Plot)))
htAll$Plot[which(is.na(htAll$Plot))] <- htAll$Plot.18[which(is.na(htAll$Plot))]

# any obviously problematic fates records
table(htAll$Survival.18,htAll$Basal.18)

table(htAll$Apical.18,htAll$Epicormic.18,htAll$Survival.18)
# 19 plants with survival = 0, and apical and/or epicormic = 1
htAll$Num[which(htAll$Survival.18==0 & (htAll$Apical.18+htAll$Epicormic.18>0))]
# changing Survival values
htAll$Survival.18[which(htAll$Survival.18==0 & (htAll$Apical.18+htAll$Epicormic.18>0))] <- 1

saveRDS(htAll,'data/hectares.rds')

# there there are lines of code that generate warnings, but they aren't errors. 
# This one is fine: In if (is.na(plot.list)) plot.list <- get.plot(year = year) :the condition has length > 1 and only the first element will be used
# the warnings about no non-missing arguments may require further examination
warnings()

