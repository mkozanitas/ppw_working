rm(list=ls())
source('scripts/11.PW_functions_local-test.R')
source('scripts/13.hectare_scripts.r')

## get hectare tree data
ht16 <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/2016/Hectares/hectare_trees_all.csv")
dim(ht16)
head(ht16)
names(ht16)
names(ht16)[which(names(ht16)=='DBH.cm')] <- 'DBH_cm'
ht16$Plot.Orig <- ht16$Plot
substr(ht16$Plot,1,3) <- 'PPW'
head(ht16)

table(ht16$Plot)
table(ht16$Year)

# it looks like trees from Plots are included in this hectare data
plot.trees <- grep('from offset',ht16$Accuracy)
length(plot.trees)
range(ht16$Num[plot.trees])
range(ht16$Num[-plot.trees],na.rm=T)

#remove from hectare data
ht16 <- ht16[-plot.trees,]
head(ht16)
dim(ht16)

#check species names
table(ht16$Species)
ht16$Species[which(ht16$Species=='QUEdec')] <- 'QUEDEC'

# now 2018 resurvey
ht18 <- read.csv('https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/2018/Hectares/Hectare_trees_all_2018.csv')
dim(ht18)
head(ht18)
names(ht18)
head(ht18$Num)
table(ht18$Plot)
ht18$Plot.Orig <- ht18$Plot
ht18$Plot <- paste('PPW',ht18$Plot,sep='')
substr(ht18$Plot,4,5) <- '13'
table(ht18$Plot)
head(ht18)
names(ht18)

# it looks like trees from original Plots are included in this hectare data
plot.trees <- which(ht18$Num<7000)
length(plot.trees)
ht18[plot.trees,]

#remove from hectare data
ht18 <- ht18[-plot.trees,]
head(ht18)
dim(ht18)

# now 2021
ht21 <- read.csv('https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/2021/Hectare_2021/Hectare_trees_all_2021.csv')
dim(ht21)
names(ht21)
ht21$Senesced <- NA
names(ht21)
ht21 <- ht21[,c(1:14,16,15)]
head(ht21)
ht21$Plot.Orig <- ht21$Plot
ht21$Plot <- 'PPW1352'
head(ht21)
str(ht21$Num)

# append 2021 data to 2018
all(names(ht18)==names(ht21))
ht18 <- rbind(ht18,ht21)
names(ht18)
str(ht18$Num)

# correct species misspellings in 2018 data
table(ht18$Species)
blSp <- which(ht18$Species=='')
ht18[blSp,c('Plot','Num')]

ht18$Species[which(ht18$Species=='QUEAGRI')] <- 'QUEAGR'
ht18$Species[which(ht18$Species=='QUEKEL ')] <- 'QUEKEL'
ht18$Species[which(ht18$Species=='QUERDOU')] <- 'QUEDOU'
ht18$Species[which(ht18$Species=='QUREGAR')] <- 'QUEGAR'
ht18$Species[which(ht18$Species=='QURELOB')] <- 'QUELOB'
ht18$Species[which(ht18$Species=='QURGAR')] <- 'QUEGAR'
table(ht18$Species)

#check for duplicate numbers
dNum <- table(ht16$Num)
dNum[which(dNum>1)]
ht16[which(ht16$Num %in% as.numeric(names(dNum[which(dNum>1)]))),c('Plot','Num')]

dNum <- table(ht18$Num)
dNum[which(dNum>1)]
tmp <- ht18[which(ht18$Num %in% as.numeric(names(dNum[which(dNum>1)]))),c('Plot','Num')]
write.csv(tmp,'data/ht18dups.csv')

# now paste together
names(ht16)
names(ht18)
names(ht16)[-3] <- paste(names(ht16)[-3],'.15',sep='')
names(ht18)[-3] <- paste(names(ht18)[-3],'.18',sep='')

dim(ht16)
dim(ht18)
str(ht16)
str(ht18)

htAll <- merge(ht16,ht18,by='Num',all = T)
dim(htAll)
names(htAll)

# species name missing in 2018
blSp <- which(htAll$Species.18=='')
length(blSp)
htAll[blSp,c('Species.15','Species.18')]

# QC - trees that have changed plots
misPlot <- (which(htAll$Plot.15!=htAll$Plot.18))
misPlot
htAll[misPlot,c('Num','Plot.15','Subplot.15','Species.15','Plot.18','Subplot.18','Species.18')]
# remove these
htAll <- htAll[-misPlot,]

# changed species
misSpp <- which(htAll$Species.15 != htAll$Species.18)
length(misSpp)
table(htAll$Species.15[misSpp],htAll$Species.18[misSpp])
tmp <- htAll[misSpp,c('Num','Plot.15','Species.15','DBH_cm.15','Plot.18','Species.18','DBH_cm.18')]
head(tmp)
write.csv(tmp,'data/ht-misSpp.csv')

# are these all senesced plants?
dim(htAll)
table(htAll$Senesced.18,useNA = 'always')


# how well do sizes match up?
plot(htAll$DBH_cm.15,htAll$DBH_cm.18,log='xy')
htAll[which.max(htAll$DBH_cm.18),]

# definitely looks like some errors, but I'll just rely on DBH.15, unless missing (mostly 2021 new superplot, presumably)
# XXXXXX flip this around??? XXXXXXX
htAll$DBH_cm.17 <- htAll$DBH_cm.15
length(which(is.na(htAll$DBH_cm.17)))
table(htAll$Plot.18[which(is.na(htAll$DBH_cm.17))])

htAll$DBH_cm.17[which(is.na(htAll$DBH_cm.17))] <- htAll$DBH_cm.18[which(is.na(htAll$DBH_cm.17))]
length(which(is.na(htAll$DBH_cm.17)))

# for now, assign 2015 species as 'correct', unless missing
htAll$Species <- htAll$Species.15
length(which(is.na(htAll$Species)))
htAll$Species[which(is.na(htAll$Species))] <- htAll$Species.18[which(is.na(htAll$Species))]

#same for plot
htAll$Plot <- htAll$Plot.15
length(which(is.na(htAll$Plot)))
htAll$Plot[which(is.na(htAll$Plot))] <- htAll$Plot.18[which(is.na(htAll$Plot))]

# now assign 2018 fates
table(htAll$Survival.18,htAll$Basal.18,useNA='always')
nosurv <- which(is.na(htAll$Survival.18))
nobasal <- which(is.na(htAll$Basal.18))
nosorb <- intersect(nosurv,nobasal)
survnobas <- which(is.na(htAll$Basal.18) & !is.na(htAll$Survival.18))

table(htAll$Plot.18[nosorb],useNA='always')
tmp <- which(!is.na(htAll$Plot.18[nosorb]))
htAll[nosorb[tmp],]

htAll[survnobas,c('Num','Plot','Survival.18','Basal.18','BasalResproutCount.18','BasalResproutHeight.18','Epicormic.18','Apical.18','Notes.18')]

table(htAll$Basal.18,htAll$BasalResproutCount.18,useNA='always')
tmp <- which(htAll$Basal.18==0 & htAll$BasalResproutCount.18==0)
length(tmp)
htAll$Num[tmp]
write.csv(htAll[tmp,],'data/basal-problems.csv')

tmp <- intersect(which(htAll$Apical.18==1 | htAll$Epicormic.18==1),which(htAll$Survival.18==0))
length(tmp)
htAll[tmp,c('Plot','Num','Survival.18','Epicormic.18','Apical.18','PercLivingCanopy.18','Senesced.18')]
write.csv(htAll[tmp,],'data/surv-apical-epi-problem.csv')

htAll$DN.18 <- NA
htAll$DN.18[which(htAll$Survival.18==0 & htAll$Basal.18==0)] <- 1
htAll$DN.18[which(htAll$Survival.18==1 | htAll$Basal.18==1)] <- 0
table(htAll$Survival.18,htAll$Basal.18,htAll$DN.18,useNA='always')

htAll$DR.18 <- NA
htAll$DR.18[which(htAll$Survival.18==0 & htAll$Basal.18==1)] <- 1
htAll$DR.18[which(htAll$Survival.18==1 | htAll$DN.18==1)] <- 0
table(htAll$Survival.18,htAll$Basal.18,htAll$DR.18,useNA='always')

htAll$LR.18 <- NA
htAll$LR.18[which(htAll$Survival.18==1 & htAll$Basal.18==1)] <- 1
htAll$LR.18[which(htAll$Survival.18==0 | htAll$Basal.18==0)] <- 0
table(htAll$Survival.18,htAll$Basal.18,htAll$LR.18,useNA='always')

htAll$LN.18 <- NA
htAll$LN.18[which(htAll$Survival.18==1 & htAll$Basal.18==0)] <- 1
htAll$LN.18[which(htAll$Survival.18==0 | htAll$LR.18==1)] <- 0
table(htAll$Survival.18,htAll$Basal.18,htAll$LN.18,useNA='always')

# check only one of the four is assigned a 1 - looks good
table(apply(htAll[,c('DN.18','DR.18','LN.18','LR.18')],1,sum),useNA='always')

htAll$fate4.18 <- NA
htAll$fate4.18[which(htAll$DN.18==1)] <- 'DN'
htAll$fate4.18[which(htAll$DR.18==1)] <- 'DR'
htAll$fate4.18[which(htAll$LN.18==1)] <- 'LN'
htAll$fate4.18[which(htAll$LR.18==1)] <- 'LR'
table(htAll$fate4.18,useNA='a')

htAll$fate3.18 <- htAll$fate4.18
htAll$fate3.18[which(htAll$fate4.18 %in% c('LN','LR'))] <- 'GC'
table(htAll$fate3.18,useNA='a')

# assign basal diameter
#D10 = DBH.cm * 1.176 + 1.070
htAll$d10.17 <- htAll$DBH_cm.17*1.176 + 1.07
plot(htAll$DBH_cm.17,htAll$d10.17)

# Assign all trees alive before the fire, all are Type='TR', and assign integer TreeNums
htAll$Live.17 <- 1
htAll$Type.17 <- 'TR'
htAll$TreeNum <- floor(htAll$Num)

# any obviously problematic fates records
table(htAll$Apical.18,htAll$Epicormic.18,htAll$Survival.18)

# 22 plants with survival = 0, and apical and/or epicormic = 1
htAll$Num[which(htAll$Survival.18==0 & (htAll$Apical.18+htAll$Epicormic.18>0))]
# changing Survival values
htAll$Survival.18[which(htAll$Survival.18==0 & (htAll$Apical.18+htAll$Epicormic.18>0))] <- 1

# 29 plants with Survival = 1 and neither Basal or Epicormic - leave for now

saveRDS(htAll,'data/hectares.rds')

