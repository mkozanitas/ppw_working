# analyze single species multinomial models
tAllm <- read.csv('data/tAll-analyzeData-update.csv',as.is=T,row.names=1)
head(tAllm)
names(tAllm)

(spAm <- sort(unique(tAllm$Species)))
(spN <- sort(table(tAllm$Species)))

(spAm <- names(spN)[which(spN>=75)])
(spAm <- spAm[-which(spAm=='BACPIL')])

# Optional - remove saplings added in 2018, where we might be introducing detection bias towards small survivors
newSap <- which(is.na(tAllm$Year.13) & tAllm$Year.18==2018 & tAllm$Type.18=='SA')
length(newSap)
tAllm$Num[newSap]
table(tAllm$Plot[newSap])

# OPTION: comment this in or out to exercise option
# tAllm <- tAllm[-newSap,]
table(tAllm$SpCd14)
dim(tAllm)

# OPTION: if needed trim data by stem size
#tAllm <- tAllm[which(tAllm$ld10>=c(-0.5)),]
dim(tAllm)


# SECTION 'THREE_FATES' - binomial and multinomial compared, for visual purposes
tAllm$Dead.18 <- 1 - tAllm$Live.18
tAllm$fate3.18 <- (-1)
tAllm$fate3.18[which(tAllm$fate.18=='DN')] <- 0
tAllm$fate3.18[which(tAllm$fate.18=='DR')] <- 1
tAllm$fate3.18[which(tAllm$fate.18 %in% c('LN','LR'))] <- 2
tAllm <- tAllm[-which(tAllm$fate3.18<0),]
table(tAllm$fate3.18) 

(f3Levels <- 0:2)
(f3PlotVals <- c(0.95,1,1.05))
f3Labs <- c('DN','DR','LR+LN')
(f3PlotCols <- c('black','red','green'))
fsPlotSize <- c(0.5,0.75,1,1.25)

tAllm$f3PlotVals <- f3PlotVals[match(tAllm$fate3.18,f3Levels)]
tAllm$f3PlotCols <- f3PlotCols[match(tAllm$fate3.18,f3PlotCols)]
tAllm$fsPlotSize <- fsPlotSize[match(tAllm$fsCat,fsPlotSize)]

# comment in or out to select one of these options
FireLevels <- c('Mod+High'); FVals <- 2:3
FireLevels <- c('Low'); FVals <- 1
FireLevels <- c('None'); FVals <- 0
FireLevels <- c('Low:High'); FVals <- 1:3 #changes FireSev range to all c(1:3)
FireLevels <- c('All'); FVals <- 0:3 #changes FireSev range to all c(1:3)

#### NEED TO SUBSET BY TYPE HERE
tAllmp <- tAllm[which(tAllm$fsCat %in% FVals & tAllm$Species %in% spAm),]
dim(tAllmp)

# remove rows with no size data
tAllmp <- tAllmp[which(!is.na(tAllmp$ld10)),]
dim(tAllmp)

# NOW FIT MULTINOMIAL
require(nnet)

# SPECIES SPECIFIC MULTINOMIAL - QUADRATIC CAN BE ADDED HERE '+ld10.2' - changes results some
spAm
spsel <- "ARBMEN"
temp <- tAllmp[which(tAllmp$Species %in% spsel),]
table(temp$fate3.18)

str(temp$fsCat)
fit.mn <- multinom(as.factor(fate3.18) ~ ld10 + ld10.2 + as.factor(fsCat), data=temp)
fit.mn 
BIC(fit.mn)
head(round(fitted(fit.mn),2))
dim(fitted(fit.mn))

plot(f3PlotVals~ld10,data=temp,col=tAllmp$f3PlotCols,ylim=c(0,1.1),main=paste('Fire Level:',FireLevels,'; Species',spsel))
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
temp$yval <- temp[,yvalname]
table(temp$yval)
plot(temp$ld10,temp$yval)

fit3 <- glm(yval~ld10+ld10.2+as.factor(fsCat),data=temp,family='binomial')
fitted(fit3)
BIC(fit3)
plot(yval~ld10,data=temp,cex=0.25,pch=19,col='blue')
for (j in 1:4) {
  rsel <- which(temp$fsCat==j-1)
  points(fitted(fit3)[rsel]~temp$ld10[rsel],cex=fsPlotSize[j])
}

