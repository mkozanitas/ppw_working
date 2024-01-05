# analyze single species multinomial models
# if anything has changed in tAll, run examineData-tAll up to this line: ## XYZABC - run to here 

rm(list=ls())
tAllm <- read.csv('data/tAll-analyzeData-update.csv',as.is=T,row.names=1)
head(tAllm)
names(tAllm)
dim(tAllm)

# remove TS
tAllm <- tAllm[which(tAllm$Type.17 %in% c('TR','SA')),]
dim(tAllm)

(spAm <- sort(unique(tAllm$Species)))
(spN <- sort(table(tAllm$Species)))

(spAm <- names(spN)[which(spN>=50)])
(spAm <- spAm[-which(spAm=='BACPIL')])

# Optional - remove saplings added in 2018, where we might be introducing detection bias towards small survivors
newSap <- which(is.na(tAllm$Year.13) & tAllm$Year.18==2018 & tAllm$Type.18=='SA')
length(newSap)
tAllm$Num[newSap]
table(tAllm$Plot[newSap])
table(tAllm$fate.18[newSap],tAllm$fsLevel[newSap])

temp <- tAllm[which(tAllm$Species %in% spAm),]
table(temp$fate.18[temp$Type.18=='SA'],temp$fsLevel[temp$Type.18=='SA'])

# find saplings in medium and high severity claimed to be LN or LR
probSA <- which(tAllm$Species %in% spAm
                & tAllm$fsLevel>1
                & tAllm$Type.18=='SA'
                & tAllm$fate.18 %in% c('LN','LR'))
length(probSA)
write.csv(tAllm[probSA,c('Plot','Num')],'data/prob-saplings.csv')

# OPTION: comment this in or out to exercise option
# tAllm <- tAllm[-newSap,]
table(tAllm$SpCd14)
dim(tAllm)

# OPTION: if needed trim data by stem size
# tAllm <- tAllm[which(tAllm$ld10>=c(0)),]
dim(tAllm)
table(tAllm$fsCat)
table(tAllm$fsLevel)

# SECTION 'THREE_FATES' - binomial and multinomial compared, for visual purposes
tAllm$Dead.18 <- 1 - tAllm$Live.18
tAllm$fate3.18 <- (-1)
tAllm$fate3.18[which(tAllm$fate.18=='DN')] <- 0
tAllm$fate3.18[which(tAllm$fate.18=='DR')] <- 1
tAllm$fate3.18[which(tAllm$fate.18 %in% c('LN','LR'))] <- 2
table(tAllm$fate3.18) 
tAllm <- tAllm[which(tAllm$fate3.18 %in% 0:2),]
table(tAllm$fate3.18) 

(f3Levels <- 0:2)
(f3PlotVals <- c(0.95,1,1.05))
f3Labs <- c('DN','DR','LR+LN')
(f3PlotCols <- c('black','red','green'))
fsPlotSize <- c(0.5,0.75,1,1.25)

tAllm$f3PlotVals <- f3PlotVals[match(tAllm$fate3.18,f3Levels)]
tAllm$f3PlotCols <- f3PlotCols[match(tAllm$fate3.18,f3Levels)]
tAllm$fsPlotSize <- fsPlotSize[match(tAllm$fsCat,fsPlotSize)]

# comment in or out to select one of these options
FireLevels <- c('Mod+High'); FVals <- 2:3
FireLevels <- c('Low'); FVals <- 1
FireLevels <- c('None'); FVals <- 0
FireLevels <- c('None:Low'); FVals <- 0:1 #changes FireSev range to all c(1:3)
FireLevels <- c('All'); FVals <- 0:3 #changes FireSev range to all c(1:3)

#### Subset to selected fire values and species
tAllmp <- tAllm[which(tAllm$fsCat %in% FVals & tAllm$Species %in% spAm),]
dim(tAllmp)

# remove rows with no size data
tAllmp <- tAllmp[which(!is.na(tAllmp$ld10)),]
dim(tAllmp)

# Translate type and size data to categories
tAllmp$SizeCat <- NA
tAllmp$SizeCat[which(tAllmp$Type.17=='SA' & tAllmp$SA.BD_cm.17 <= max(tAllmp$SA.BD_cm.17,na.rm=T))] <- 'SA2'
#tAllmp$SizeCat[which(tAllmp$Type.17=='SA' & tAllmp$SA.BD_cm.17<1)] <- 'SA1'
tAllmp$SizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 <= max(tAllmp$DBH_cm.17,na.rm=T))] <- 'TR4'
#tAllmp$SizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 < 40)] <- 'TR3'
tAllmp$SizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 < 20)] <- 'TR2'
tAllmp$SizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 < 10)] <- 'TR1'
table(tAllmp$SizeCat,useNA='always')

head(which(is.na(tAllmp$SizeCat)))
tAllmp[563,]

# NOW FIT MULTINOMIAL
require(nnet)

# SPECIES SPECIFIC MULTINOMIAL - QUADRATIC CAN BE ADDED HERE '+ld10.2' - changes results some
spAm
TSp <- c('QUEKEL','QUEDOU','ARBMEN','QUEGAR','QUEAGR','PSEMEN','UMBCAL')
HSp <- c('QUEKEL','QUEDOU','ARBMEN','QUEGAR','QUEAGR','UMBCAL','QUEBER')
WOSp <- c('QUEDOU','QUEGAR')
RSHSp <- c('HETARB','AMOCAL','FRACAL')
NRSHSp <- 'ARCMAN'

# CHOOSE ONE
spsel <- spAm;spname <- 'All'
spsel <- TSp;spname <- 'Trees'
spsel <- HSp;spname <- 'Hardwoods'
spsel <- RSHSp;spname <- 'Resprouting Shrubs'
spsel <- NRSHSp;spname <- 'Non-Resprouting Shrub (ARCMAN)'
spsel <- WOSp;spname <- 'White Oaks'
spsel <- "UMBCAL";spname <- spsel

temp <- tAllmp[which(tAllmp$Species %in% spsel),]
nrow(temp)
table(temp$fate3.18)
table(temp$SizeCat,useNA='always')

# multinomial discrete size
temp <- temp[which(!is.na(temp$SizeCat)),]
dim(temp)
table(temp$SizeCat)

(tot <- table(temp$SizeCat,temp$fsCat))
barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,main=paste(spname,'Mortality'))
barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,main=paste(spname,'Topkill-Resprout'))
barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,main=paste(spname,'Green-Crown'))

(tot <- table(temp$fsCat,temp$SizeCat))
barplot(table(temp$fsCat,temp$SizeCat,temp$fate3.18)[,,1]/tot,beside = T,main=paste(spname,'Mortality'))
barplot(table(temp$fsCat,temp$SizeCat,temp$fate3.18)[,,2]/tot,beside = T,main=paste(spname,'Topkill-Resprout'))
barplot(table(temp$fsCat,temp$SizeCat,temp$fate3.18)[,,3]/tot,beside = T,main=paste(spname,'Green-Crown'))


fit.mn <- multinom(as.factor(fate3.18) ~ as.factor(SizeCat) + as.factor(fsCat), data=temp)
fit.mn 
BIC(fit.mn)
head(round(fitted(fit.mn),2))
dim(fitted(fit.mn))
plot(fitted(fit.mn)[,3]~ld10,data=temp)



(tot <- table(temp$SizeCat,temp$fsCat))
round(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1],2)
round(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,2)
round(tapply(fitted(fit.mn)[,1],list(temp$SizeCat,temp$fsCat),mean,na.rm=T),2)

(tot <- table(temp$SizeCat,temp$fsCat))
round(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2],2)
round(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,2)
round(tapply(fitted(fit.mn)[,2],list(temp$SizeCat,temp$fsCat),mean,na.rm=T),2)

(tot <- table(temp$SizeCat,temp$fsCat))
round(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3],2)
round(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,2)
round(tapply(fitted(fit.mn)[,3],list(temp$SizeCat,temp$fsCat),mean,na.rm=T),2)

plot(f3PlotVals~ld10,data=temp,col=temp$f3PlotCols,ylim=c(0,1.1),main=paste('Fire Level:',FireLevels,'; Species',spsel),pch=19)
i=1
j=1
for (i in 1:3) for (j in 1:4) {
  rsel <- which(temp$fsCat==j-1)
  points(temp$ld10[rsel],fitted(fit.mn)[rsel,i],col=f3PlotCols[i],cex=fsPlotSize[j])
}


# multinomial continuous size
fit.mn <- multinom(as.factor(fate3.18) ~ poly(ld10,2):as.factor(fsCat), data=temp)
fit.mn 
BIC(fit.mn)
head(round(fitted(fit.mn),2))
dim(fitted(fit.mn))

plot(f3PlotVals~ld10,data=temp,col=temp$f3PlotCols,ylim=c(0,1.1),main=paste('Fire Level:',FireLevels,'; Species',spsel),pch=19)
i=1
j=1
for (i in 1:3) for (j in 1:4) {
  rsel <- which(temp$fsCat==j-1)
  points(temp$ld10[rsel],fitted(fit.mn)[rsel,i],col=f3PlotCols[i],cex=fsPlotSize[j])
}



#binomial
#set yvalue
# For some outcomes, like NR, need to drop some species
yvalname <- 'gCrown.18'
if (yvalname=='gCrown.18') pcol <- 'green' else if (yvalname=='DR.18') pcol <- 'red' else pcol='black'
temp$yval <- temp[,yvalname]
table(temp$yval)
#plot(temp$ld10,temp$yval)

fit3 <- glm(yval~poly(ld10,2) + as.factor(fsCat),data=temp,family='binomial')
summary(fit3)
fitted(fit3)
BIC(fit3)
plot(yval~ld10,data=temp)
for (j in 1:4) {
  rsel <- which(temp$fsCat==j-1)
  points(fitted(fit3)[rsel]~temp$ld10[rsel],cex=fsPlotSize[j],col=pcol,pch=19)
}