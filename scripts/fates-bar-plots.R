# make bar plots by species and functional group
# if anything has changed in tAll, run examineData-tAll up to this line: ## XYZABC - run to here 

rm(list=ls())
tAllm <- read.csv('data/tAll-analyzeData-update.csv',as.is=T,row.names=1)
head(tAllm)
names(tAllm)
dim(tAllm)

# remove TS
tAllm <- tAllm[which(tAllm$Type.17 %in% c('TR','SA')),]
dim(tAllm)

# read in species attribute table
spAtt <- read.csv('data/all-spp-attributes.csv',row.names = 1)
head(spAtt)
nrow(spAtt)

# extract list of species in data
(spA <- sort(unique(tAllm$Species)))

# check all species in data are in list
all(spA %in% spAtt$Species)

# make match var to extract attributes from table
s2s <- match(spA,spAtt$Species)

(spN <- sort(table(tAllm$Species)))

# make list of common species for plotting
(spAc <- names(spN)[which(spN>=50)])
(spAc <- spAc[-which(spAc %in% c('BACPIL','QUEBER'))])
(spAc <- sort(spAc))

# Optional - remove saplings added in 2018, where we might be introducing detection bias towards small survivors
newSap <- which(is.na(tAllm$Year.13) & tAllm$Year.18==2018 & tAllm$Type.18=='SA')
length(newSap)
tAllm$Num[newSap]
table(tAllm$Plot[newSap])
table(tAllm$fate.18[newSap],tAllm$fsLevel[newSap])

temp <- tAllm[which(tAllm$Species %in% spAc),]
table(temp$fate.18[temp$Type.18=='SA'],temp$fsLevel[temp$Type.18=='SA'])

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
tAllmp <- tAllm[which(tAllm$fsCat %in% FVals & tAllm$Species %in% spA),]
dim(tAllmp)

# remove rows with no size data
tAllmp <- tAllmp[which(!is.na(tAllmp$ld10)),]
dim(tAllmp)

# Translate type and size data to categories
tAllmp$TSizeCat <- NA
tAllmp$TSizeCat[which(tAllmp$Type.17=='SA' & tAllmp$SA.BD_cm.17 <= max(tAllmp$SA.BD_cm.17,na.rm=T))] <- 'SA'
#tAllmp$SizeCat[which(tAllmp$Type.17=='SA' & tAllmp$SA.BD_cm.17<1)] <- 'SA1'
tAllmp$TSizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 <= max(tAllmp$DBH_cm.17,na.rm=T))] <- 'TR3'
tAllmp$TSizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 < 20)] <- 'TR2'
tAllmp$TSizeCat[which(tAllmp$Type.17=='TR' & tAllmp$DBH_cm.17 < 10)] <- 'TR1'
table(tAllmp$TSizeCat,useNA='always')

tAllmp$SSizeCat <- NA
tAllmp$SSizeCat[which(tAllmp$Type.17=='SA')] <- 'SA'
tAllmp$SSizeCat[which(tAllmp$Type.17=='TR')] <- 'TR'
table(tAllmp$SSizeCat,useNA='always')

# Plot for all trees and hardwoods species
tsets <- list()
tset.names <- c()
tsets[[1]] <- spA[which(spAtt$Shrub.Tree[s2s]=="T")]
tset.names[1] <- 'All trees'
tsets[[2]] <- tsets[[1]][-which(tsets[[1]]=='PSEMEN')]
tset.names[2] <- 'Hardwood trees'
tsets[[3]] <- intersect(tsets[[1]],spAc)
tset.names[3] <- 'Common trees'
tsets[[4]] <- intersect(tsets[[2]],spAc)
tset.names[4] <- 'Common hardwoods'

i=1
for (i in 1:length(tsets))
{
  spsel <- tsets[[i]]
  spname <- tset.names[i]
  temp <- tAllmp[which(tAllmp$Species %in% spsel),]
  nrow(temp)
  
  temp <- temp[which(temp$TSizeCat %in% c('SA','TR1','TR2','TR3')),]
  temp$SizeCat <- temp$TSizeCat
  nrow(temp)
  
  (tot <- table(temp$SizeCat,temp$fsCat))
  
  pdf(paste('fates-figs/fates-',spname,'.pdf',sep=''),6,6)
  op=par(mfrow=c(2,2))
  
  tree.cols <- c('grey90','grey60','grey30','black')
  barplot(tot,beside=T,col=tree.cols,main=paste(spname,'Sample Sizes'))
  barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
  barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=tree.cols,main=paste(spname,'Topkill-Resprout'),ylim=c(0,1))
  barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,col=tree.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
  
  par(op)
  dev.off()
}

# Plot for all shrubs and resprouting species
ssets <- list()
sset.names <- c()
ssets[[1]] <- spA[which(spAtt$Shrub.Tree[s2s]=="S")]
ssets[[1]] <- ssets[[1]][-which(ssets[[1]]=='BACPIL')]
sset.names[1] <- 'All shrubs'
ssets[[2]] <- ssets[[1]][-which(ssets[[1]]=='ARCMAN')]
sset.names[2] <- 'Resprouting shrubs'
ssets[[3]] <- intersect(ssets[[1]],spAc)
sset.names[3] <- 'Common shrubs'
ssets[[4]] <- intersect(ssets[[2]],spAc)
sset.names[4] <- 'Common resprouting shrubs'

i=1
for (i in 1:4)
{
  spsel <- ssets[[i]]
  spname <- sset.names[i]
  temp <- tAllmp[which(tAllmp$Species %in% spsel),]
  nrow(temp)
  
  temp <- temp[which(temp$SSizeCat %in% c('SA','TR')),]
  temp$SizeCat <- temp$SSizeCat
  nrow(temp)    
  
  (tot <- table(temp$SizeCat,temp$fsCat))
  
  pdf(paste('fates-figs/fates-',spname,'.pdf',sep=''),6,6)
  op=par(mfrow=c(2,2))
  
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  barplot(tot,beside=T,col=shrub.cols,main=paste(spname,'Sample Sizes'))
  barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=shrub.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
  barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=shrub.cols,main=paste(spname,'Topkill-Resprout'),ylim=c(0,1))
  barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,col=shrub.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
  
  par(op)
  dev.off()
}

spsel <- ssets[[2]]
spname <- sset.names[2]
temp <- tAllmp[which(tAllmp$Species %in% spsel),]
# one of the surprising results here is high green crown for shrubs in high severity fire. How many were there?
table(temp$Species[which(temp$fsCat==3)]) #13

# which ones had green crown?
table(temp$Species[which(temp$fsCat==3 & temp$fate3.18=='2')]) # 1 CEACUN, 4 CEAPAR, 2 SAMNIG

# should these be recoded as topkill-resprout? Or might some (CEACUN) have germinated after the fire. Here are the numbers and plots
# is it possible these regrew in summer following fire, rather than surviving. Certainly seems likely for CEAPAR. If so, they should be removed, like BACPIL.
temp[which(temp$fsCat==3 & temp$fate3.18=='2'),c('Plot','Species')]

i <- 6
# Start loop here, through individual species
for (i in 1:length(spAc)) 
{
  spsel <- spAc[i];spname <- spsel
  temp <- tAllmp[which(tAllmp$Species %in% spsel),]
  type <- spAtt$Shrub.Tree[which(spsel==spAtt$Species)]
  nrow(temp)
  if (type=='S')
  {
    temp <- temp[which(temp$SSizeCat %in% c('SA','TR')),]
    temp$SizeCat <- temp$SSizeCat
    nrow(temp)    
  } else
  {
    temp <- temp[which(temp$TSizeCat %in% c('SA','TR1','TR2','TR3')),]
    temp$SizeCat <- temp$TSizeCat
    nrow(temp)
  }
  (tot <- table(temp$SizeCat,temp$fsCat))
  
  pdf(paste('fates-figs/fates-',spsel,'.pdf',sep=''),6,6)
  op=par(mfrow=c(2,2))
  if (type=='S') 
  {
    shrub.cols <- c('grey90','black')
    barplot(tot,beside=T,col=shrub.cols,main=paste(spname,'Sample Sizes'))
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=shrub.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
    if (spsel=='ARCMAN') {
      plot(1:10,1:10,xlab='',ylab='',type='n',xaxt='n',yaxt='n')
      text(5,5,'Non-sprouting species',cex=1)
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=shrub.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
    } else {
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=shrub.cols,main=paste(spname,'Topkill-Resprout'),ylim=c(0,1))
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,col=shrub.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
    }
  } else 
  {
    tree.cols <- c('grey90','grey60','grey30','black')
    barplot(tot,beside=T,col=tree.cols,main=paste(spname,'Sample Sizes'))
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
    if (spsel=='PSEMEN') {
      plot(1:10,1:10,xlab='',ylab='',type='n',xaxt='n',yaxt='n')
      text(5,5,'Non-sprouting species',cex=1)
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=tree.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
    } else {
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=tree.cols,main=paste(spname,'Topkill-Resprout'),ylim=c(0,1))
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,col=tree.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
    }
  }
  par(op)
  dev.off()
}



## now run multinomial models for hardwoods and resprouting shrubs. What terms are justified to include
require(nnet)

# Common Hardwood Trees
(spsel <- tsets[[4]])
(spname <- tset.names[4])
temp <- tAllmp[which(tAllmp$Species %in% spsel),]
nrow(temp)

temp$fsCat <- as.factor(temp$fsCat)
temp$TSizeCat <- as.factor(temp$TSizeCat)

fit.mn0 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat + Species + northness, data=temp)
BIC(fit.mn0) #4049

fit.mn1 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat + Species, data=temp)
BIC(fit.mn1) # 4062, clearly fit.mn0 is better, keep northness

fit.mn2 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat + northness, data=temp)
BIC(fit.mn2) # 4013 Much better, drop species

fit.mn3 <- multinom(as.factor(fate3.18) ~ TSizeCat + Species + northness, data=temp)
BIC(fit.mn3) # 5551 fsCat has huge effect

fit.mn4 <- multinom(as.factor(fate3.18) ~ fsCat + Species + northness, data=temp)
BIC(fit.mn4) # 4876 TSizeCat has huge effect

# so, start with model 2
fit.mn2 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat + northness, data=temp)
BIC(fit.mn2) # 4014

# add interactions one at a time
fit.mn.i1 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat  + northness + TSizeCat:fsCat, data=temp)
BIC(fit.mn.i1) # 4029 No, don't keep

fit.mn.i2 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat  + northness + TSizeCat:northness, data=temp)
BIC(fit.mn.i2) # 4014 Equal, don't need

fit.mn.i3 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat  + northness + fsCat:northness + TSizeCat:fsCat, data=temp)
BIC(fit.mn.i3) # 4009 Yes, keep

fit.mn.i4 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat  + northness + fsCat:northness + TSizeCat:fsCat + fsCat:northness:fsCat, data=temp)
BIC(fit.mn.i4) # 4009 Same, third order not needed

# Interaction model 3 supported
fit.mn.i3 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat  + northness + fsCat:northness + TSizeCat:fsCat, data=temp)
BIC(fit.mn.i3) # 4009
coefficients(fit.mn.i3)
write.csv(coefficients(fit.mn.i3),'fates-figs/coeff-common-hardwoods.csv')

# Model for PSEMEN
(spsel <- "PSEMEN")
(spname <- "PSEMEN")
temp <- tAllmp[which(tAllmp$Species %in% spsel),]
nrow(temp)

temp$fsCat <- as.factor(temp$fsCat)
temp$TSizeCat <- as.factor(temp$TSizeCat)

fit.mn0 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat + northness, data=temp)
BIC(fit.mn0) #698

fit.mn1 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat, data=temp)
BIC(fit.mn1) # 690, better, drop northness

fit.mn3 <- multinom(as.factor(fate3.18) ~ TSizeCat + northness, data=temp)
BIC(fit.mn3) # 1034 fsCat has huge effect

fit.mn4 <- multinom(as.factor(fate3.18) ~ fsCat + northness, data=temp)
BIC(fit.mn4) # 972 TSizeCat has huge effect

# so, start with model 1
fit.mn1 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat, data=temp)
BIC(fit.mn1) # 691, better, drop northness

# add interactions one at a time
fit.mn.i1 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat  + TSizeCat:fsCat, data=temp)
BIC(fit.mn.i1) # 739 No, don't keep

# No interaction models are strongly supported as better, so go with simplest: model 1
fit.mn1 <- multinom(as.factor(fate3.18) ~ TSizeCat + fsCat, data=temp)
BIC(fit.mn1) # 699, better, drop northness
coefficients(fit.mn1)
write.csv(t(coefficients(fit.mn1)),'fates-figs/coeff-psemen.csv')

# Common Resprouting Shrubs
(spsel <- ssets[[4]])
(spname <- sset.names[4])
temp <- tAllmp[which(tAllmp$Species %in% spsel),]
nrow(temp)

temp$fsCat <- as.factor(temp$fsCat)
temp$TSizeCat <- as.factor(temp$TSizeCat)

fit.mn0 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat + Species + northness, data=temp)
BIC(fit.mn0) #2008

fit.mn1 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat + Species, data=temp)
BIC(fit.mn1) # 1994, Lower, drop northness

fit.mn2 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat + northness, data=temp)
BIC(fit.mn2) # 2064 Much worse, keep species

fit.mn3 <- multinom(as.factor(fate3.18) ~ SSizeCat + Species + northness, data=temp)
BIC(fit.mn3) # 2369 fsCat has huge effect

fit.mn4 <- multinom(as.factor(fate3.18) ~ fsCat + Species + northness, data=temp)
BIC(fit.mn4) # 2000 Not supported relative to model 1

fit.mn5 <- multinom(as.factor(fate3.18) ~ fsCat + Species, data=temp)
BIC(fit.mn5) #1987 - best model

# so, start with model 5

# add interaction
fit.mn.i1 <- multinom(as.factor(fate3.18) ~ fsCat*Species, data=temp)
BIC(fit.mn.i1) # 2033 No, don't keep


# No interaction models are strongly supported as better, so go with simplest: model 5
fit.mn5 <- multinom(as.factor(fate3.18) ~ fsCat + Species, data=temp)
BIC(fit.mn5) #1987 - best model
coefficients(fit.mn5)
write.csv(coefficients(fit.mn5),'fates-figs/coeff-resp-shrubs.csv')

# Model for ARCMAN
(spsel <- "ARCMAN")
(spname <- "ARCMAN")
temp <- tAllmp[which(tAllmp$Species %in% spsel),]
nrow(temp)
table(temp$fsCat)

# one plant in high severity, not enough to test effect
temp <- temp[-which(temp$fsCat==3),]
temp$fsCat <- as.numeric(temp$fsCat)
temp$fsCat <- as.factor(temp$fsCat)
str(temp$fsCat)
temp$SSizeCat <- as.factor(temp$SSizeCat)

fit.mn0 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat + northness, data=temp)
BIC(fit.mn0) #70.4

fit.mn1 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat, data=temp)
BIC(fit.mn1) # 66.5, better, drop northness

fit.mn3 <- multinom(as.factor(fate3.18) ~ SSizeCat + northness, data=temp)
BIC(fit.mn3) # 81.96 fsCat has huge effect

fit.mn4 <- multinom(as.factor(fate3.18) ~ fsCat + northness, data=temp)
BIC(fit.mn4) # 66.4 SSizeCat has huge effect

fit.mn5 <- multinom(as.factor(fate3.18) ~ fsCat , data=temp)
BIC(fit.mn5) # 62.7 Lowest BIC - just fire severity

fit.mn6 <- multinom(as.factor(fate3.18) ~ SSizeCat , data=temp)
BIC(fit.mn6) # 82.7 Lowest BIC - just fire severity

coefficients(fit.mn5)

# so, start with model 2
fit.mn1 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat, data=temp)
BIC(fit.mn1) # 67, better, drop northness

# add interactions one at a time
fit.mn.i1 <- multinom(as.factor(fate3.18) ~ SSizeCat + fsCat  + SSizeCat:fsCat, data=temp)
BIC(fit.mn.i1) # 73.2 No, don't keep

# No interaction models are strongly supported as better, so go with simplest: model 2
fit.mn5 <- multinom(as.factor(fate3.18) ~ fsCat , data=temp)
BIC(fit.mn5) # 62.7 Lowest BIC - just fire severity
coefficients(fit.mn5)
write.csv(t(coefficients(fit.mn5)),'fates-figs/coeff-ARCMAN.csv')

table(temp$fsCat)
