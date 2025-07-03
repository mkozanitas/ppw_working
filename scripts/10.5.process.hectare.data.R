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
table(ht16$Subplot)

# examine in context +/- 1
pm1 <- function(x) c(x-1,x,x+1)
cb <- which(ht16$Subplot=='C')
cb1 <- 854:865
ht16[cb1,]

table(ht16$Plot)
table(ht16$Year)

# any duplicate numbers:
dn <- table(ht16$Num)
dn[which(dn>1)]

# it looks like trees from Plots are included in this hectare data
plot.trees <- grep('from offset',ht16$Accuracy)
length(plot.trees)
range(ht16$Num[plot.trees])
range(ht16$Num[-plot.trees],na.rm=T)

#remove from hectare data
dim(ht16)
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
ht18$Plot <- paste('PPW',ht18$Plot,sep='')
ht18$Plot.Orig <- ht18$Plot
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
ht21$Plot.Orig <- paste('PPW',ht21$Plot,sep='')
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
dNum[which(dNum>1)] # NONE !!!!!
#ht16[which(ht16$Num %in% as.numeric(names(dNum[which(dNum>1)]))),c('Plot','Num')]

dNum <- table(ht18$Num)
dNum[which(dNum>1)]
tmp <- ht18[which(ht18$Num %in% as.numeric(names(dNum[which(dNum>1)]))),c('Plot','Num')]
tmp
#write.csv(tmp,'data/ht18dups.csv')

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

# species name missing in 2018  -NONE
# blSp <- which(htAll$Species.18=='')
# length(blSp) # NONE
# htAll[blSp,c('Species.15','Species.18')]

# QC - trees that have changed plots
misPlot <- (which(htAll$Plot.15!=htAll$Plot.18))
misPlot # NONE
htAll[misPlot,c('Num','Plot.15','Subplot.15','Species.15','Plot.18','Subplot.18','Species.18')]
# remove these
htAll <- htAll[-misPlot,]

# changed species - NONE!
misSpp <- which(htAll$Species.15 != htAll$Species.18)
length(misSpp)
# table(htAll$Species.15[misSpp],htAll$Species.18[misSpp])
# tmp <- htAll[misSpp,c('Num','Plot.15','Species.15','DBH_cm.15','Plot.18','Species.18','DBH_cm.18')]
# tmp
# write.csv(tmp,'data/ht-misSpp.csv')

# are these all senesced plants?
# dim(htAll)
# table(htAll$Senesced.18,useNA = 'always')

# how well do sizes match up?
plot(htAll$DBH_cm.15,htAll$DBH_cm.18,log='xy')
htAll[which.max(htAll$DBH_cm.18),]

# definitely looks like some errors, but I'll just rely on DBH.15, unless missing (mostly 2021 new superplot, presumably)

# some 2018 had 0 dbh instead of NA
tmp0 <- which(htAll$DBH_cm.18==0)
htAll$DBH_cm.18[tmp0] <- NA

htAll$DBH_cm.17 <- htAll$DBH_cm.18
length(which(is.na(htAll$DBH_cm.17)))
table(htAll$Plot.18[which(is.na(htAll$DBH_cm.17))])

htAll$DBH_cm.17[which(is.na(htAll$DBH_cm.17))] <- htAll$DBH_cm.15[which(is.na(htAll$DBH_cm.17))]
length(which(is.na(htAll$DBH_cm.17)))

# for now, assign 2015 species as 'correct', unless missing
tmp <- which(is.na(htAll$Species.15) & htAll$Plot.18!='PPW1352')
length(tmp)
htAll[tmp,c('Num','Plot.15','Species.15','Plot.18','Species.18','DBH_cm.18')]
htAll$Species <- htAll$Species.15
length(which(is.na(htAll$Species)))
htAll$Species[which(is.na(htAll$Species))] <- htAll$Species.18[which(is.na(htAll$Species))]
length(which(is.na(htAll$Species)))

#same for plot
htAll$Plot <- htAll$Plot.15
length(which(is.na(htAll$Plot)))
htAll$Plot[which(is.na(htAll$Plot))] <- htAll$Plot.18[which(is.na(htAll$Plot))]

# how many missing in 2015
a15p18 <- which(is.na(htAll$Plot.15) & !is.na(htAll$Plot.18) & htAll$Plot.18!='PPW1352')
length(a15p18)
htAll[a15p18,c('Num','Plot','Species','DBH_cm.17','Plot.18','Species.18','DBH_cm.18')]
write.csv(htAll[a15p18,],'data/abs15pres18.csv')

# how many missing in 2018
tmp <- which(!is.na(htAll$Plot.15) & is.na(htAll$Plot.18))
length(tmp)
htAll[tmp,c('Num','Plot.15','Species.15','DBH_cm.15','Plot.18','Species.18','DBH_cm.18')]
write.csv(htAll[tmp,],'data/pres15abs18.csv')

# how many under 20 cm
dbhlt20 <- which(htAll$DBH_cm.17<20)
length(dbhlt20)
range(htAll$DBH_cm.17[dbhlt20],useNA=T)
hist(htAll$DBH_cm.17[dbhlt20])

# check species again and remove species not included in hectares
# Hetarb are too small anyway, so don't exclude anything else
table(htAll$Species[which(htAll$DBH_cm.17>=20)])

#Create ExcStem flag for analysis step
htAll$ExcStem <- 0
htAll$ExcStem[which(htAll$Species=='HETARB')] <- 1
htAll$ExcStem[dbhlt20] <- 1

# now assign 2018 fates
table(htAll$Survival.18,htAll$Basal.18,useNA='always')
nosurv <- which(is.na(htAll$Survival.18))
nobasal <- which(is.na(htAll$Basal.18))
nosorb <- intersect(nosurv,nobasal)
survnobas <- which(is.na(htAll$Basal.18) & !is.na(htAll$Survival.18))
htAll[survnobas,]
htAll[survnobas,'ExcStem'] <- 1

table(htAll$Plot.18[nosorb],useNA='always')
tmp <- which(!is.na(htAll$Plot.18[nosorb]))
htAll[nosorb[tmp],]
htAll[nosorb[tmp],'ExcStem'] <- 1

# tmp <- intersect(which(htAll$Apical.18==1 | htAll$Epicormic.18==1),which(htAll$Survival.18==0))
# length(tmp)
# htAll[tmp,c('Plot','Num','Survival.18','Epicormic.18','Apical.18','PercLivingCanopy.18','Senesced.18')]
# write.csv(htAll[tmp,],'data/surv-apical-epi-problem.csv')

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

# who is missing
tmp <- which(is.na(htAll$fate3.18))
htAll[tmp,c('Num','Plot','Species','DBH_cm.17','Survival.18','Basal.18','ExcStem')]

# fate of trees missing in 2015 and found in 2018
htAll[a15p18,c('Num','Species','Plot','DBH_cm.17','fate3.18','ExcStem')]

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

tmp <- which(htAll$Survival.18==1 & (htAll$Apical.18+htAll$Epicormic.18==0))
htAll[tmp,c('Num','Plot.18','PercLivingCanopy.18')]
length(tmp)
write.csv(htAll[tmp,],'data/Surv1ApEpPer0.csv')

# Assign subplots based on coordinates
noUTM15 <- which(is.na(htAll$UTM.E.15) | is.na(htAll$UTM.N.15))
length(noUTM15)

noUTM18 <- which(is.na(htAll$UTM.E.18) | is.na(htAll$UTM.N.18))
length(noUTM18)

noUTM <- intersect(noUTM15,noUTM18)
table(htAll$ExcStem[noUTM])

htAll$UTM.E.17 <- htAll$UTM.E.15
htAll$UTM.E.17[which(is.na(htAll$UTM.E.17))] <- htAll$UTM.E.18[which(is.na(htAll$UTM.E.17))]
htAll$UTM.N.17 <- htAll$UTM.N.15
htAll$UTM.N.17[which(is.na(htAll$UTM.N.17))] <- htAll$UTM.N.18[which(is.na(htAll$UTM.N.17))]

plot(htAll$UTM.E.17,htAll$UTM.N.17,asp=1)

# read plot coordinates
htc <- read.csv('input_data/plot_info/hectares-18-20m-FS.csv')
head(htc)
htc$Plot.Orig <- htc$Plot
substr(htc$Plot,4,5) <- '13'

htAll$Subplot.17 <- 'XX'
htAll$coordOutlier <- 0

drawBox <- function(b,exlim=50,plot=NULL) {
  plot(b[1],b[3],xlim=c(b[1]-exlim,b[2]+exlim),ylim=c(b[3]-exlim,b[4]+exlim),pch=19,cex=2,col='red',main=plot)
  points(b[1],b[4],pch=19,cex=2,col='red')
  points(b[2],b[3],pch=19,cex=2,col='red')
  points(b[2],b[4],pch=19,cex=2,col='red')
}

plot.list <- sort(unique(htc$Plot))
i=1
for (i in 1:length(plot.list)) {
  plot <- plot.list[i]
  rs <- which(htAll$Plot==plot)
  tmp <- htAll[rs,c('Plot','Subplot.17','Num','UTM.E.17','UTM.N.17')]
  plot(htAll$UTM.E.17,htAll$UTM.N.17,asp=1,type='n',main=tmp[1,'Plot'])
  points(tmp$UTM.E.17,tmp$UTM.N.17,pch=19)
  
  cs <- which(htc$Plot==plot)
  uxm <- min(htc$UTM.x[cs])
  uxx <- max(htc$UTM.x[cs]+20)
  uym <- min(htc$UTM.y[cs])
  uyx <- max(htc$UTM.y[cs]+20)
  b <- c(uxm,uxx,uym,uyx)
  
  #drawBox(b,exlim=50,plot)
  #points(tmp[,c('UTM.E.17','UTM.N.17')],asp=1)
  source('scripts/adjustHectCoords.R')
  #drawBox(b,exlim=50,plot)
  #points(tmp$UTM.E.17+as.numeric(optRes[1]),tmp$UTM.N.17+as.numeric(optRes[2]))
  print(optRes)
  
  tmp$UTM.E.17 <- tmp$UTM.E.17 + as.numeric(optRes[1])
  tmp$UTM.N.17 <- tmp$UTM.N.17 + as.numeric(optRes[2])
  
  nrow(tmp)
  tmp$Subplot.17 <- 'L0'
  j=103
  for (j in 1:nrow(tmp)) {
    if (is.na(tmp$UTM.E.17[j] | tmp$UTM.N.17[j])) tmp$Subplot.17[j] <- NA else {
      if (tmp$UTM.E.17[j] >= htc$UTM.x[cs[1]]) substr(tmp$Subplot.17[j],1,1) <- 'A'
      if (tmp$UTM.E.17[j] >= htc$UTM.x[cs[6]]) substr(tmp$Subplot.17[j],1,1) <- 'B'
      if (tmp$UTM.E.17[j] >= htc$UTM.x[cs[11]]) substr(tmp$Subplot.17[j],1,1) <- 'C'
      if (tmp$UTM.E.17[j] >= htc$UTM.x[cs[16]]) substr(tmp$Subplot.17[j],1,1) <- 'D'
      if (tmp$UTM.E.17[j] >= htc$UTM.x[cs[21]]) substr(tmp$Subplot.17[j],1,1) <- 'E'
      if (tmp$UTM.E.17[j] > htc$UTM.x[cs[21]]+20) substr(tmp$Subplot.17[j],1,1) <- 'M'
      if (tmp$UTM.N.17[j] >= htc$UTM.y[cs[1]]) substr(tmp$Subplot.17[j],2,2) <- '1'
      if (tmp$UTM.N.17[j] >= htc$UTM.y[cs[2]]) substr(tmp$Subplot.17[j],2,2) <- '2'
      if (tmp$UTM.N.17[j] >= htc$UTM.y[cs[3]]) substr(tmp$Subplot.17[j],2,2) <- '3'
      if (tmp$UTM.N.17[j] >= htc$UTM.y[cs[4]]) substr(tmp$Subplot.17[j],2,2) <- '4'
      if (tmp$UTM.N.17[j] >= htc$UTM.y[cs[5]]) substr(tmp$Subplot.17[j],2,2) <- '5'
      if (tmp$UTM.N.17[j] > htc$UTM.y[cs[5]]+20) substr(tmp$Subplot.17[j],2,2) <- '9'
    }
  }
  table(tmp$Subplot.17)
  htAll[rs,c('Plot','Subplot.17','Num','UTM.E.17','UTM.N.17')] <- tmp
}
table(htAll$Subplot.17)
htAll$coordOutlier[which(htAll$Subplot.17 %in% c('A0','A9','B0','B9','C0','C9','D0','D9','E0','E9') | substr(htAll$Subplot.17,1,1) %in% c('L','M'))] <- 1
table(htAll$coordOutlier,htAll$ExcStem)

# hardcode long outliers as Exclude for Analysis, pending Melina corrections
htAll$ExcStem[which(htAll$Num %in% c(10568,10573,10998,99200,11046.1,11083,11083.1,99191,11158,11158.1,11158.2,99220))] <- 1

# Now move near outliers to nearest plot
table(htAll$Subplot.17)
htAll$Subplot.17[which(htAll$Subplot.17=='A0')] <- 'A1'
htAll$Subplot.17[which(htAll$Subplot.17=='B0')] <- 'B1'
htAll$Subplot.17[which(htAll$Subplot.17=='C0')] <- 'C1'
htAll$Subplot.17[which(htAll$Subplot.17=='D0')] <- 'D1'
htAll$Subplot.17[which(htAll$Subplot.17=='E0')] <- 'E1'

htAll$Subplot.17[which(htAll$Subplot.17=='A9')] <- 'A5'
htAll$Subplot.17[which(htAll$Subplot.17=='B9')] <- 'B5'
htAll$Subplot.17[which(htAll$Subplot.17=='C9')] <- 'C5'
htAll$Subplot.17[which(htAll$Subplot.17=='D9')] <- 'D5'
htAll$Subplot.17[which(htAll$Subplot.17=='E9')] <- 'E5'

htAll$Subplot.17[which(htAll$Subplot.17=='L1')] <- 'A1'
htAll$Subplot.17[which(htAll$Subplot.17=='L2')] <- 'A2'
htAll$Subplot.17[which(htAll$Subplot.17=='L3')] <- 'A3'
htAll$Subplot.17[which(htAll$Subplot.17=='L4')] <- 'A4'
htAll$Subplot.17[which(htAll$Subplot.17=='L5')] <- 'A5'

htAll$Subplot.17[which(htAll$Subplot.17=='M1')] <- 'E1'
htAll$Subplot.17[which(htAll$Subplot.17=='M2')] <- 'E2'
htAll$Subplot.17[which(htAll$Subplot.17=='M3')] <- 'E3'
htAll$Subplot.17[which(htAll$Subplot.17=='M4')] <- 'E4'
htAll$Subplot.17[which(htAll$Subplot.17=='M5')] <- 'E5'

htAll$Subplot.17[which(htAll$Subplot.17=='L0')] <- 'A1'
htAll$Subplot.17[which(htAll$Subplot.17=='L9')] <- 'E1'
htAll$Subplot.17[which(htAll$Subplot.17=='M0')] <- 'A5'
htAll$Subplot.17[which(htAll$Subplot.17=='M9')] <- 'E5'
table(htAll$Subplot.17,htAll$ExcStem)

plot(htAll$UTM.E.17,htAll$UTM.N.17,asp=1)
points(htAll[which(htAll$coordOutlier==1),c('UTM.E.17','UTM.N.17')],col='red')

# identifyPlot# identify this datanoUTM18# identify this data set as hectares
htAll$Survey <- 'Hect'
head(htAll)

saveRDS(htAll,'data/hectares.rds')

