#analyzeData-tAll is to compare 2013, 2018, 2019, and 2020 for delayed mortality and other changes in state

# RUN prepareData and examineData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R

rm(list=ls())
library(lme4)
library(glmmTMB)

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

spNames <- read.csv('data/all-spp-names.csv',row.names = 1)
head(spNames)
names(spNames)[which(names(spNames)=='x')] <- 'spName'

# input merged dataframe
tAll <- read.csv('data/tAll.csv',as.is=T)

# initial growth analysis to identify potentially problematic size data
# suggest taking median eliminating most negative and very large outliers. Suggests median diameter growth of 0.2 cm in five years
tAll$ddbh.1318 <- tAll$dbh.18-tAll$dbh.13
hist(tAll$ddbh.1318[tAll$UseForBAGrowth])
summary(tAll$ddbh.1318[tAll$UseForBAGrowth])

plot(tAll$dbh.13,tAll$ddbh.1318)
abline(h=0)

tAll$absddbh.1318 <- abs(tAll$ddbh.1318)
hist(tAll$absddbh.1318)
tail(sort(tAll$absddbh.1318))

write.csv(tAll[which(tAll$absddbh.1318>5),c('Num','Plot.18')],'data/growthquestions.csv')

## END CREATE PROXY prefire data for new trees encountered post-fire, including new plots

## THIS CODE SECTION ANALYZES 2018 FATES ###

# ANALYZE BY SPECIES AND TYPE for 2018 post-fire fates
(use.species <- spNames$spName)

# create fst dataframe - FateSummaryTable
fst12 <- data.frame(Species=rep(use.species,each=2),Type=rep(c('SA','TR'),length(use.species)),N13=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
head(fst12)                  
tail(fst12)

i=5
for (i in 1:nrow(fst12))
{
  sp <- fst12$Species[i]
  ty <- fst12$Type[i]
  temp <- tAll[which(tAll$Species.13==sp & tAll$Type.13==ty),]
  
  fst12$N13[i] <- sum(temp$Live.13,na.rm=T)
  
  ## The next three lines are all equivalent - just using third one
  #fst12$N18.DN[i] <- length(which(temp$fate.18=='DN'))
  #fst12$N18.DN[i] <- length(which(temp$DN.18=='1'))
  fst12$N18.DN[i] <- sum(temp$DN.18,na.rm = T)
  
  fst12$N18.DR[i] <- sum(temp$DR.18,na.rm = T)
  fst12$N18.LN[i] <- sum(temp$LN.18,na.rm = T)
  fst12$N18.LR[i] <- sum(temp$LR.18,na.rm = T)
  miss <- which(temp$Live.13==1 & is.na(temp$DN.18)==1)
  if (length(miss)>0) for (j in 1:length(miss)) print(temp[miss[j],c('Plot.13','Num')])
  fst12$nMissing <- fst12$N13-(fst12$N18.DN+fst12$N18.DR+fst12$N18.LN+fst12$N18.LR)
}

fst12
head(fst12)
tail(fst12)
sum(fst12$nMissing) #fixed 1330 dups we were previously ignoring- should now be zero

# Summary across all species and types #GHYTIGYFYGIGH - all categories at once in a table 
tree.sum <- apply(fst12[fst12$Type=='TR',-c(1:2)],2,sum)
(tree.sum-tree.sum[6])/(tree.sum[1]-tree.sum[6])

sap.sum <- apply(fst12[fst12$Type=='SA',-c(1:2)],2,sum)
(sap.sum-sap.sum[6])/(sap.sum[1]-sap.sum[6])

# NaN errors?
#ts.sum <- apply(fst12[fst12$Type=='TS',-c(1:2)],2,sum,na.rm=T)
#(ts.sum-sap.sum[6])/(ts.sum[1]-ts.sum[6])

fate.sum <- apply(fst12[,-c(1:2)],2,sum)
(fate.sum-fate.sum[6])/(fate.sum[1]-fate.sum[6])

#GHYTIGYFYGIGH - same as above but just for one category at a time (DN/DR/LR/LN)

SArows <- which(fst12$Type=='SA')
TRrows <- which(fst12$Type=='TR')
#TSrows <- which(fst12$Type=='TS') 

fst12$percSurv <- 1 - fst12$N18.DN/fst12$N13
fst12$percSurv[fst12$N13==0] <- NA

#copied AbSp code from below to calculate percSurv of each species in each category

## now run by species for species with lots of data
spN <- table(tAll$Species.13)
spN[order(spN)]

(spA <- names(spN)[which(spN>=25)])
spA <- spA[-which(spA=='QUEBER')] # drop QUEBER - none died!

# dropping QUEBER, FRACAL, AMOCAL for now
AbSp <- c('ARBMEN','ARCMAN','HETARB','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')
length(AbSp)

# use abundant species identified above for ESA
(spA <- AbSp)

fst12a <- fst12[which(fst12$Species %in% AbSp),]

# look at TR and SA in order of percent Survival
fst12a[fst12a$Type=='TR',][order(fst12a$Species[fst12a$Type=='TR']),]
fst12a[fst12a$Type=='SA',][order(fst12a$Species[fst12a$Type=='SA']),]
#fst12a[fst12a$Type=='TS',][order(fst12a$percSurv[fst12a$Type=='TS']),]


## choose fire severity metric
fsmet <- 'Tubbs.MTBS.RDNBR.30'
names(fs)
summary(fs[,fsmet])
hist(fs[,fsmet])
sort(fs[,fsmet])

f2t <- match(tAll$Plot.13,fs$Plot)
head(f2t)
tail(f2t)
tAll$FireSev <- fs[f2t,fsmet]
dim(tAll)
tail(tAll)
table(tAll$Plot.13)

# RDNBR fire severity levels
# Unburned
# Low <170, but burned
# Medium 170-700
# High: >700
# unburned: 1308, 1309, 1311, 1312, 1327, 1344, 1347

# Create a discrete FireSev variable
summary(tAll$FireSev)
hist(tAll$FireSev,breaks=c(-100,0,50,100,150,200,300,400,500,600,700,800,1000))
fs.breaks <- c(-100,135,430,1000) # Based on Parks et al. 2014, but not using intermediate split at 304
tAll$fsCat <- cut(tAll$FireSev,fs.breaks)
table(tAll$fsCat)
tAll$fsLevel <- as.numeric(tAll$fsCat)
table(tAll$fsLevel)

# manually code unburned as 0
unburned.plots <- c('PPW1308','PPW1309','PPW1311','PPW1312','PPW1327','PPW1344','PPW1347')
tAll$fsLevel[which(tAll$Plot.13 %in% unburned.plots)] <- 0
table(tAll$fsLevel)
tAll$fsCat <- as.factor(tAll$fsLevel)
table(tAll$fsCat)
str(tAll$fsCat)

# plots experiencing each fire severity level, and how many N13 individuals in each
table(tAll$Plot.13[which(tAll$fsLevel==0)])
table(tAll$Plot.13[which(tAll$fsLevel==1)])
table(tAll$Plot.13[which(tAll$fsLevel==2)])
table(tAll$Plot.13[which(tAll$fsLevel==3)])

# # Initial examination of dbh
# dim(tAll)
# length(which(tAll$Type.18=='TS'))
# nodbh <- which(is.na(tAll$dbh.13))
# length(nodbh)
# head(tAll[nodbh,])
# summary(tAll$dbh.13)

## Use d10 for analysis - remember for SA is original basal data and for TR is calculated from dbh
# D10 = DBH.cm * 1.176 + 1.070

# if want to change and use size from a different year, change here. Then from here on ld10 is generic
hist(tAll$d10.18)
tAll$ld10 <- log10(tAll$d10.13)
tAll$ld10[which(!is.finite(tAll$ld10))] <- NA
hist(tAll$ld10)
summary(tAll$ld10,useNA='always')

# now move diameters forward from 2013, if it's missing in 2018 
tAll$ld10[which(is.na(tAll$ld10))] <- log10(tAll$d10.18[which(is.na(tAll$ld10))])
summary(tAll$ld10,useNA='always')
hist(tAll$ld10)

# smallest tree now had d10 = 1*1.176+1.07. log10 of this = 0.3514

# create ts to switch between types
types <- c('SA','TR')
op=par(mfrow=c(1,2))
for (i in 1:2) {
  ty <- types[i]
  print(hist(tAll$ld10[tAll$Type.18==ty],main=paste('Type',ty)))
}
par(op)

# create squared variable for quadratic analysis
tAll$ld10.2 <- tAll$ld10^2

# check on fates
table(tAll$Live.18, tAll$fate.18)

# now try saplings and trees separately
op <- par(mfrow=c(1,2))
for (i in 1:2) plot(Live.18~ld10,data=tAll[tAll$Type.18==types[i],],main=paste('Type',types[i]))
par(op)

# or plot together
plot(Live.18~ld10,data=tAll,main=paste('All',types[i]))

# check fate values
table(tAll$Resprout.18,tAll$fate.18)

#### SURVIVAL ANALYSIS - first cut, size only!!
types <- c('TR','SA')

# adjust values here to subset data, for TR and/or SA and fire severity level and species
tAlls <- tAll[which(tAll$Type.18 %in% types[1:2] & tAll$fsLevel>=0 & tAll$Species.18 %in% AbSp),]

# if needed trim data by stem size
#tAlls <- tAlls[which(tAlls$ld10>=c(-0.5)),]

# sample size
tAlls$fPlot <- as.factor(tAlls$Plot.18)
nrow(tAlls)
(N <- length(which(!is.na(tAlls$ld10) & !is.na(tAlls$Live.18))))

#set yvalue
yvalname <- 'LN.18'
tAlls$yval <- tAlls[,yvalname]

#fit model
fit0 <- glm(yval~ld10,data=tAlls,family='binomial')
BIC(fit0)

fit1 <- glm(yval~ld10+ld10.2,data=tAlls,family='binomial')
#summary(fit1)
AIC(fit1)
BIC(fit1)

fit1n <- glm(yval~ld10+ld10.2+northness,data=tAlls,family='binomial')
#summary(fit1)
AIC(fit1n)
BIC(fit1n)
coefficients(fit1n)

# we can check all combinations, but for Live.18 and the full data, northness wasn't justified to add to the full model, with FS and species. So for now we're including it, but need to check, for each variable, before reporting in the paper

fit2 <- glmer(yval~ld10+ld10.2+northness+(1|fPlot),data=tAlls,family='binomial')
#summary(fit2)
AIC(fit2)
BIC(fit2)
fitPlots <- rownames(coefficients(fit2)$fPlot)

# now add fire severity
# check whether model can fit fsLevel and random plot factor - converges!
fit3 <- glm(yval~ld10+ld10.2+northness+fsCat,data=tAlls,family='binomial')
#summary(fit3)
AIC(fit3)
BIC(fit3)

fit4 <- glmer(yval~ld10+ld10.2+northness+fsCat+(1|fPlot),data=tAlls,family='binomial')
#summary(fit4)
BIC(fit4)
# YES! (for both types, absp, all fsLevels....)

fit5 <- glm(yval~ld10+ld10.2+fsCat+northness+Species.18,data=tAlls,family='binomial')
BIC(fit5)
coefficients(fit5)

# here's the full model with FS and species, and no northness - check via BIC whether it can be included in final results
fit5x <- glm(yval~ld10+ld10.2+fsCat+Species.18,data=tAlls,family='binomial')
BIC(fit5x)

fit5l <- glm(yval~ld10+fsCat+northness+Species.18,data=tAlls,family='binomial')
BIC(fit5l)

#fit6 <- glmer(yval~ld10+ld10.2+fsCat+Species.18+(1|fPlot),data=tAlls,family='binomial')
# DOESN'T CONVERGE combining species and random factor plots

# made newdata for prediction
nvals <- 11
ld10vals <- seq(min(tAlls$ld10,na.rm=T),max(tAlls$ld10,na.rm=T),length.out=nvals)
ld10vals <- c(min(tAlls$ld10,na.rm=T),log10(c(1,1.5,2,5,10,15,20,50,100)),max(tAlls$ld10,na.rm=T))

nd <- expand.grid(ld10vals,fitPlots,unique(tAlls$fsCat),AbSp)
names(nd) <- c('ld10','fPlot','fsCat','Species.18')
nd$fPlot <- as.factor(nd$fPlot)
nd$fsCat <- as.factor(nd$fsCat)
nd$ld10.2 <- nd$ld10^2
nd$northness <- 0
head(nd)
dim(nd)

#xx <- data.frame(ld10=seq(min(tAlls$ld10,na.rm=T),max(tAlls$ld10,na.rm=T),length.out=nvals))
#nd2 <- data.frame(fPlot=rep(as.factor(fitPlots),nrow(xx)),ld10=rep(xx$ld10,each=length(fitPlots)))
#nd3 <- data.frame(fPlot=rep(as.factor(fitPlots),nrow(nd2)),ld10=rep(nd2$ld10,nrow(nd2)),fsCat <- )
#rm('xx')


# predict value from fit
nd$predVal0 <- predict(fit0,newdata=nd,type='response')
nd$predVal1 <- predict(fit1,newdata=nd,type='response')
nd$predVal2 <- predict(fit2,newdata=nd,type='response')
nd$predVal3 <- predict(fit3,newdata=nd,type='response')
nd$predVal4 <- predict(fit4,newdata=nd,type='response')
nd$predVal5 <- predict(fit5,newdata=nd,type='response')
nd$predVal5l <- predict(fit5l,newdata=nd,type='response')
head(nd)
#plot(nd$predVal3,nd$predVal4)
#plot(nd$predVal2,nd$predVal4)
#plot(nd$predVal4,nd$predVal5)

#plot data and predicted values
range(tAlls$ld10,na.rm=T)
plot(tAlls$ld10,tAlls$yval,main=yvalname)
points(nd$ld10,nd$predVal0,lwd=4) # linear size
points(nd$ld10,nd$predVal1,lwd=4) # quadratic size
points(nd$ld10,nd$predVal2,lwd=4) # qsize + plots only
points(nd$ld10,nd$predVal3,lwd=4) # qsize + fire levels
points(nd$ld10,nd$predVal4,lwd=4) # qsize + fire levels + plots
points(nd$ld10,nd$predVal5,lwd=4) # qsize + fire levels + species
points(nd$ld10,nd$predVal5l,lwd=4) # qsize + fire levels + species

#rowSel <- which(nd$Species.18=='UMBCAL' & nd$fsCat==1)
#points(nd$ld10[rowSel],nd$predVal5[rowSel],lwd=4,col='red') # qsize + fire levels + species

pAvg <- rep(NA,11)
i=1
fs=1
for (fsv in 1:4) {
  fs <- fsv-1
  for (i in 1:nvals) {
    ld <- ld10vals[i]
    pAvg[i] <- mean(nd$predVal5[which(nd$ld10==ld & nd$fsCat==fs)])
  }
  #print(pAvg)
  points(ld10vals,pAvg,col='red',lwd=4,type='b')
}
nd[which(nd$ld10==0&nd$fsCat==0),]

#abline(h=0.5,lty=2) 
#h draws horizontal line at .5, lty -dashed or solid, lwd is line width

## We have an issue with apparent higher survival in small individuals. Let's investigate
#plot(Live.18~ld10,data=tAlls)

# next section not edited to separate TR and SA, not valid!!
if (FALSE) {
  #what is critical basal area to achieve 50% survival?
  (ld50 <- nd$ld10[which(nd$pSurvAll>=0.5)[1]])
  abline(v=ld50,lty=2)
  
  (spRes <- data.frame(species='All',N=N,ld50=ld50,slp=fit$coefficients[2]))
  spRes

  
  i <- 3
  for (i in 1:length(spA))
  {
    (species <- spA[i])
    tAlls <- tAll[which(tAll$Species.13==species),]
    N <- length(which(!is.na(tAlls$Basal.Area.13) & !is.na(tAlls$Live.18)))
    fit <- glm(Live.18~ld10,data=tAlls,family='binomial')
    pval <- predict(fit,newdata=nd,type='response')
    nd <- data.frame(nd,pval)
    names(nd)[length(names(nd))] <- paste('pSurv_',species,sep='')
    lines(nd$ld10,nd[,ncol(nd)])
    
    ld50 <- NA
    if (fit$coefficients[2]>0) 
      if (nd[1,ncol(nd)]<0.5) 
        ld50 <- nd$ld10[which(nd[,ncol(nd)]>=0.5)[1]] 
    spRes <- rbind(spRes,c(species,N,ld50,fit$coefficients[2]))
  }
  spRes$ld50 <- as.numeric(spRes$ld50)
  spRes$d50 <- round(10^spRes$ld50,3)
  spRes$ld50 <- round(as.numeric(spRes$ld50),3)
  spRes$slp <- round(as.numeric(spRes$slp),3)
  spRes
  
  plotSP <- function(tAll,species=NULL,xlims=range(tAll$ld10,na.rm=T))
  {
    tAlls <- tAll[which(tAll$Species.13==species),]
    mindbh <- min(tAlls$ld10,na.rm=T)
    maxdbh <- max(tAlls$ld10,na.rm=T)
    plot(tAlls$ld10,tAlls$Live.18,main=species,xlim=xlims)
    fit <- glm(Live.18~ld10,data=tAlls,family='binomial')
    nd <- with(tAll,data.frame(ld10=seq(mindbh,maxdbh,length.out=1001)))
    pval <- predict(fit,newdata=nd,type='response')
    lines(nd$ld10,pval)  
  }
  plotSP(tAll,'ARCMAN')
  plotSP(tAll,'QUEAGR')
  plotSP(tAll,'QUEDOU')
  plotSP(tAll,'ARBMEN')
  plotSP(tAll,'HETARB')
  plotSP(tAll,'PSEMEN')
  plotSP(tAll,'QUEGAR')
  plotSP(tAll,'QUEKEL')
  plotSP(tAll,'UMBCAL')
}

########### Starting with Models########

## Full model with species (can subset out data with mod/high fslevel(2:3) for visualization only- FS not in model-change back to c(0:3 to expand to all fslevels) to generate curves

# remember variable types
types
ty <- types[1:2]
fsl <- 0:3
tAlls <- tAll[which(tAll$Type.18 %in% ty & tAll$Species.18 %in% spA & tAll$fsLevel %in% fsl),]
dim(tAlls)

# assign dependent variable to rVar # can switch between Live.18 and gCrown.18 to generate figures
names(tAlls)
selVar <- 'gCrown.18'
tAlls$rVar <- tAlls[,selVar]

(N <- length(which(!is.na(tAlls$ld10) & !is.na(tAlls$rVar))))
table(tAlls$Species.13)
#head(tAlls)
fit1 <- glm(rVar~ld10 + ld10.2 + Species.18,data=tAlls,family='binomial')
fit2 <- glmer(rVar~ld10 + ld10.2 + Species.18 + (1|fPlot),data=tAlls,family='binomial')
# fails to converge

table(tAlls$Species.18,tAlls$Plot.18)
summary(fit)

summary(tAlls$ld10)
#choose sizes here at which predicted values will be calculated, can change size 1 to 2cm and rerun
predSizes <- c(10,30,80)
nd <- with(tAlls,data.frame(ld10=rep(log10(predSizes),length(spA)),Species.13=rep(spA,each=length(predSizes))))
nd$ld10.2 <- nd$ld10^2
nd
nd$pValue <- round(predict(fit,newdata=nd,type='response'),4)
#- can use Sp.names instead of Species.13- below to reorder for barplot
nd$Sp.Names <- c("aARCMAN","bPSEMEN", "cQUEDOU", "dQUEKEL", "eARBMEN", "fQUEGAR", "gUMBCAL", "hHETARB", "iQUEAGR")


# with shrubs
selSize <- 3 #change size 1 above from 10cm to 2cm in order to look at shrubs
barplot(pValue~Species.13,data=nd[which(nd$ld10==log10(predSizes[selSize])),],ylim=c(0,1),main=paste('Type',ty,'FSLev',fsl,selVar,'predval@',predSizes[selSize],'cm dbh'))

# remove shrubs for larger sizes cm plot
# op=par will stack 3 figures (1,3) for horizontal & (3,1) for vertical
op=par(mfrow=c(1,3)) 
selSize <- 1
barplot(pValue~Species.13,data=nd[which(nd$ld10==log10(predSizes[selSize])  & !nd$Species.13 %in% c('HETARB','ARCMAN')),],ylim=c(0,1),main=paste('Type',ty,'FSLev',fsl,selVar,'predval@',predSizes[selSize],'cm dbh'))

selSize <- 2
barplot(pValue~Species.13,data=nd[which(nd$ld10==log10(predSizes[selSize])  & !nd$Species.13 %in% c('HETARB','ARCMAN')),],ylim=c(0,1),main=paste('Type',ty,'FSLev',fsl,selVar,'predval@',predSizes[selSize],'cm dbh'))

selSize <- 3
barplot(pValue~Species.13,data=nd[which(nd$ld10==log10(predSizes[selSize])  & !nd$Species.13 %in% c('HETARB','ARCMAN')),],ylim=c(0,1),main=paste('Type',ty,'FSLev',fsl,selVar,'predval@',predSizes[selSize],'cm dbh'))


#run to reset to single figure 
par(op)


### MODELS WITH FIRE SEVERITY - use 4 levels as factors
# Change response variable here and then run model. Swap gCrown.18 or Live.18 to analyze crown survival
ty <- 'TR'
tAlls <- tAll[which(tAll$Type.18==ty & tAll$Species.13 %in% spA & tAll$fsLevel %in% c(0:3)),]

#did a fix here to set variable like we did above using selVar..
selVar <- 'DN.18'
tAlls$rVar <- tAlls[,selVar]
fit <- glm(rVar~ld10 + as.factor(fsLevel),data=tAlls,family='binomial')
summary(fit)

nvals <- 101
ld10.vals <- seq(min(tAlls$ld10,na.rm=T),max(tAlls$ld10,na.rm=T),length.out=nvals)
#FireSev.vals <- seq(min(fs[,fsmet],na.rm=T),max(fs[,fsmet],na.rm=T),length.out=nvals)
FireLevel.vals <- sort(unique(tAlls$fsLevel))
nd <- with(tAll,data.frame(ld10=rep(ld10.vals,length(FireLevel.vals)),fsLevel=rep(FireLevel.vals,each=nvals)))
head(nd)

nd$rVarPred <- predict(fit,newdata=nd,type='response')
head(nd)
tail(nd)

# Plot response variable as a function of size, with isoclines as a function fire severity. Lines are increasing - survival is higher for larger trees, but lower at higher fire severity
plot(tAlls$ld10,tAlls$rVar,type="n", xlab="log10DBH", ylab = selVar)
i=1
for (i in 1:length(FireLevel.vals)) {
  ndt <- nd[which(nd$fsLevel==FireLevel.vals[i]),]
  lines(ndt$ld10,ndt$rVarPred,col="red")
  text(ld10.vals[50],ndt$rVarPred[which(ndt$ld10==ld10.vals[50])],FireLevel.vals[i])
}

  
#### THIS WON'T WORK NOW AS MODEL ABOVE WAS CHANGED TO USE FIRE SEVERITY LEVELS
# Plot survival as a function of fire severity, with isoclines as a function ld10. Lines are declining - survival is lower at higher fire severity, but larger trees have higher values
# plot(tAlls$FireSev,tAlls$Live.18)
# i=1
# for (i in 1:nvals) {
#   ndt <- nd[which(nd$ld10==ld10.vals[i]),]
#   lines(ndt$FireSev,ndt$pSurvAll)
#   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ld10.vals[i])
# }
### END COMMENT OUT

## RECODED WITH FIRE LEVELS
# model with size, fire, species, predicting survival
{
  fit2 <- glm(rVar~ld10 + as.factor(fsLevel) + Species.13,data=tAlls,family='binomial')
  fit1 <- glm(rVar~ld10 + as.factor(fsLevel) + Species.13+ as.factor(fsLevel):Species.13,data=tAlls,family='binomial')
  BIC(fit1)
  BIC(fit2)
  summary(fit1)
  summary(fit2)
  
  nvals <- 11
  ld10.vals <- seq(min(tAlls$ld10,na.rm=T),max(tAlls$ld10,na.rm=T),length.out=nvals)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ld10=rep(ld10.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ld10=rep(nd$ld10,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  plotSpecies <- function(spname,nd.tAlls=nd2) {
    tAlls <- tAlls[which(tAlls$Species.13==spname),]
    #op=par(mfrow=c(1,2))
    
    # plot(tAlls$FireSev,tAlls$Live.18,main=spname,xlim=range(tAlls$FireSev,na.rm=T))
    # i=1
    # for (i in 1:nvals) {
    #   ndt <- nd2[which(nd2$Species.13==spname & nd2$ld10==ld10.vals[i]),]
    #   lines(ndt$FireSev,ndt$pSurvAll)
    #   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ld10.vals[i])
    # }
    
    plot(tAlls$ld10,tAlls$Live.18,xlim=range(tAlls$ld10,na.rm=T))
    i=2
    for (i in 1:length(unique(nd.tAlls$fsLevel))) {
      ndt <- nd.tAlls[which(nd2$Species.13==spname & nd.tAlls$fsLevel==fsLevels.vals[i]),]
      lines(ndt$ld10,ndt$pSurvAll)
      text(ld10.vals[5],ndt$pSurvAll[which(ndt$ld10==ld10.vals[5])],fsLevels.vals[i])
    }
    #par(op)
  }
  #This will plot isoclines of survival at each fire severity for any indv species 
  plotSpecies('PSEMEN')
  plotSpecies('ARBMEN')
  plotSpecies('QUEAGR')
  plotSpecies('UMBCAL')
  plotSpecies('ARBMEN')
  
  # 4/26/23 - everything working to here!!!!
  
  #using data from model w/ fire sev included- not subsetted out for visualization as above
  # obtain predicted main effects of fire severity for each species at a common size
  nvals <- 1
  ld10.vals <- log10(2)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ld10=rep(ld10.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ld10=rep(nd$ld10,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  
  ###these plots don't illustrate survival in mod+high fire severity as well as previous model (excluding FS) with data subsetted out for FS levels
  
  #barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==0),],ylim=c(0,1))
  #barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==1),],ylim=c(0,1))
  #barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==2),],ylim=c(0,1))
  #barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==3),],ylim=c(0,1))
}
# same code, with gCrown.18 instead of Live.18
{
  fit2 <- glm(gCrown.18~ld10 + as.factor(fsLevel) + Species.13,data=tAlls,family='binomial')
  fit1 <- glm(gCrown.18~ld10 + as.factor(fsLevel) + Species.13+ as.factor(fsLevel):Species.13,data=tAlls,family='binomial')
  BIC(fit1)
  BIC(fit2)
  summary(fit1)
  summary(fit2)
  
  nvals <- 11
  ld10.vals <- seq(-0.5,2,length.out=nvals)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ld10=rep(ld10.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ld10=rep(nd$ld10,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  spname='QUEAGR'
  plotSpecies <- function(spname,nd.tAlls=nd2) {
    tAlls <- tAlls[which(tAlls$Species.13==spname),]
    #op=par(mfrow=c(1,2))
    
    # plot(tAlls$FireSev,tAlls$Live.18,main=spname,xlim=range(tAlls$FireSev,na.rm=T))
    # i=1
    # for (i in 1:nvals) {
    #   ndt <- nd2[which(nd2$Species.13==spname & nd2$ld10==ld10.vals[i]),]
    #   lines(ndt$FireSev,ndt$pSurvAll)
    #   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ld10.vals[i])
    # }
    
    plot(tAlls$ld10,tAlls$gCrown.18,xlim=range(tAlls$ld10,na.rm=T))
    i=2
    for (i in 1:length(unique(nd.tAlls$fsLevel))) {
      ndt <- nd.tAlls[which(nd2$Species.13==spname & nd.tAlls$fsLevel==fsLevels.vals[i]),]
      lines(ndt$ld10,ndt$pSurvAll)
      text(ld10.vals[5],ndt$pSurvAll[which(ndt$ld10==ld10.vals[5])],fsLevels.vals[i])
    }
    #par(op)
  }
  plotSpecies('QUEGAR')
  
  # obtain predicted main effects of fire severity for each species at a common size (2CM)
  nvals <- 1
  ld10.vals <- log10(2)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ld10=rep(ld10.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ld10=rep(nd$ld10,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  ###these plots don't illustrate survival in mod+high fire severity as well as previous model (excluding FS) with data subsetted out for FS levels 
  # barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==0),],ylim=c(0,1))
  # barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==1),],ylim=c(0,.5))
  # barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==2),],ylim=c(0,.25))
  # barplot(pSurvAll ~ Species.13,data=nd2[which(nd2$fsLevel==3),],ylim=c(0,.05))
}
## END RECODE HERE ##

#changed TB.18 to DR.18 from here on... i think thats right bc TB meant topkilled with basal 

tAll$Dead.18 <- 1 - tAll$Live.18
names(tAll)
table(tAll$Dead.18)
table(tAll$Live.18)
table(tAll$DR.18)
table(tAll$gCrown.18)
table(tAll$Live.18,tAll$DR.18)
table(tAll$Live.18,tAll$gCrown.18)

tAll$PFstatus.18 <- (-1)
tAll$PFstatus.18[which(tAll$Live.18==0)] <- 0
tAll$PFstatus.18[which(tAll$DR.18==1)] <- 1
tAll$PFstatus.18[which(tAll$gCrown.18==1)] <- 2
table(tAll$PFstatus.18) 
#this table totals 6945 out of total 6946 trees, why is one indv missing??


(PFstatusLevels <- 0:2)
(PFsPlotVals <- c(0.95,1,1.05))
(PFsPlotCols <- c('black','red','green'))

tAll$PFsPlotVals <- PFsPlotVals[match(tAll$PFstatus,PFstatusLevels)]
tAll$PFsPlotCols <- PFsPlotCols[match(tAll$PFstatus,PFstatusLevels)]

# range of ld10 for abundant species
spArows <- which(tAll$Species.13 %in% spA)
(tAllld10.range <- c(min(tAll$ld10[spArows],na.rm=T),max(tAll$ld10[spArows],na.rm=T)))

# head(tAll$fsLevel)
# head(match(tAll$fsLevel,PFstatusLevels))
# head(tAll$PFsPlotVals)
## ANALYSIS FOR ONE SPECIES
spA

# pick one of these!
selSpecies <- spA # use spA for all abundant species, rather than one species
# comment in or out next two lines
FireLevels <- c('Mod+High'); FVals <- 2:3
#FireLevels <- c('ANY LEVEL'); FVals <- 1:3 #changes FireSev range to all c(1:3)

#### NEED TO SUBSET BY TYPE HERE
tAllsp <- tAll[which(tAll$Species.13 %in% c(selSpecies) & tAll$fsLevel %in% FVals),] #individual species?
{
  #tAllsp <- tAll[which(tAll$Species.13 %in% spA),] #abundant species?
  #tAllsp <- tAll[which(tAll$Species.13 %in% spA & tAll$fsLevel>1),] #abundant sp with fs level of 1 or more?
  
  dim(tAllsp)
  tAllsp <- tAllsp[which(!is.na(tAllsp$ld10)),]
  dim(tAllsp)
  
  # code to run a binomial model and plot response curve with data
  # MORTALITY
  nd <- with(tAllsp,data.frame(ld10=seq(min(tAllsp$ld10,na.rm=T),max(tAllsp$ld10,na.rm=T),length.out=101)))
  nd$ld10.2 <- nd$ld10^2
  
  #plot(tAllsp$ld10,tAllsp$Dead.18,xlim=tAllld10.range)
  fit1 <- glm(Dead.18~ld10,data=tAllsp,family='binomial')
  fit2 <- glm(Dead.18~ld10+ld10.2,data=tAllsp,family='binomial')
  BIC(fit1)
  BIC(fit2)
  nd$pMortality <- predict(fit1,newdata=nd,type='response')
  #lines(nd$ld10,nd$pMortality)
  
  #TOPKILL WITH RESPROUT
  #plot(tAllsp$ld10,tAllsp$TB,xlim=tAllld10.range)
  fit <- glm(DR.18~ld10+ld10.2 ,data=tAllsp,family='binomial')
  summary(fit)
  nd$pDR <- predict(fit,newdata=nd,type='response')
  #lines(nd$ld10,nd$pTB)
  
  #GREEN CROWN
  #plot(tAllsp$ld10,tAllsp$gCrown.18,xlim=tAllld10.range)
  fit1 <- glm(gCrown.18~ld10,data=tAllsp,family='binomial')
  fit2 <- glm(gCrown.18~ld10+ld10.2,data=tAllsp,family='binomial')
  summary(fit1)
  nd$pGCrown <- predict(fit1,newdata=nd,type='response')
  #lines(nd$ld10,nd$pGCrown)
  
  #plot all three (change main from selSpecies to ".." to alter main title)
  plot(tAllsp$ld10,tAllsp$PFsPlotVals,col=tAllsp$PFsPlotCols,pch=19,ylim=c(-0.05,1.05),xlim=tAllld10.range,main=paste("ALL_SP",FireLevels))
  points(tAllsp$ld10,rep(-0.05,length(tAllsp$ld10)))
  lines(nd$ld10,nd$pGCrown,col='green')
  lines(nd$ld10,nd$pDR,col='red')
  lines(nd$ld10,nd$pMortality)
  
  # does the sum of the three binomials for these three exclusive fates sum to 1?
  nd$pTOT <- apply(nd[,c('pGCrown','pDR','pMortality')],1,sum)
  summary(nd$pTOT)
}

#reset plotting window with one panel
par(mfrow=c(1,1))

# NOW FIT MULTINOMIAL
require(nnet)

# MULTINOMIAL - QUADRATIC CAN BE ADDED HERE '+ld10.2' - changes results some
fit1 <- multinom(as.factor(PFstatus.18) ~ ld10 +ld10.2, data=tAllsp)
fit1
head(round(fitted(fit1),2))
dim(fitted(fit1))

plot(tAllsp$ld10,tAllsp$Live.18)
points(tAllsp$ld10,fitted(fit1)[,1],col='black')
points(tAllsp$ld10,fitted(fit1)[,2],col='red')
points(tAllsp$ld10,fitted(fit1)[,3],col='green')
summary(apply(fitted(fit1)[,1:3],1,sum))
#######################

# END CLEAN CODE HERE!!! redundant stuff happening below


###old figures- this doesn't work anymore###

# spA
# png('figures/logit1.png',width = 800, height = 1200)
# op=par(mfrow=c(4,2))
# plotSpecies('AMOCAL')
# plotSpecies('ARBMEN')
# plotSpecies('ARCMAN')
# plotSpecies('FRACAL')
# par(op)
# dev.off()
# 
# png('figures/logit2.png',width = 800, height = 1200)
# op=par(mfrow=c(4,2))
# plotSpecies('HETARB')
# plotSpecies('PSEMEN')
# plotSpecies('UMBCAL')
# par(op)
# dev.off()
# 
# png('figures/logit3.png',width = 800, height = 1200)
# op=par(mfrow=c(4,2))
# plotSpecies('QUEAGR')
# plotSpecies('QUEDOU')
# plotSpecies('QUEGAR')
# plotSpecies('QUEKEL')
# par(op)
# dev.off()

### full model with random plot
#fit <- glmer(Live.18~ld10 + FireSev + Species.13 + (1 |Plot.13),data=tAlls,family='binomial')

######## TOPKILL ANALYSIS
# current scoring has dead trees as topkilled. Change so topkill is only for those that are alive
tAll$TopkillLive.18 <- 0
tAll$TopkillLive.18[which(tAll$Topkill.18==1 & tAll$Live.18==1)] <- 1

table(tAll$Dead.18,tAll$TopkillLive.18,tAll$gCrown.18)
tAll[which(tAll$gCrown.18==0 & tAll$TopkillLive.18==0 & tAll$Dead.18==0),]
# TWO INDIVIDUALS WITH PROBLEM DATA: 3415, 4437 (already identified those above - if they've been fixed and don't show up at this point, delete this line)

tAlla <- tAll[which(tAll$FireSev>100),]

dim(tAlla)
names(tAlla)
head(tAlla)
tAlla$ld10.2 <- tAlla$ld10^2

op=par(mfrow=c(1,1))
plot(Dead.18~ld10,data=tAlla)
fitD <- glm(Dead.18~ld10+ld10.2,data=tAlla,family='binomial')
fitD
nd <- data.frame(ld10=seq(min(tAlla$ld10,na.rm=T),max(tAlla$ld10,na.rm=T),length.out=101))
nd$ld10.2 <- nd$ld10^2
head(nd)
nd$pDead <- predict(fitD,nd,type='response')
lines(pDead~ld10,data=nd,lwd=2,col='black')

#plot(Topkill.18~ld10,data=tAlla)
fitT <- glm(TopkillLive.18~ld10+ld10.2,data=tAlla,family='binomial')
fitT
nd$pTopKill <- predict(fitT,nd,type='response')
lines(pTopKill~ld10,data=nd,lwd=2,col='red')

#plot(gCrown.18~ld10,data=tAlla)
fitG2 <- glm(gCrown.18~ld10+ld10.2,data=tAlla,family='binomial')
fitG2
BIC(fitG2)
fitG1 <- glm(gCrown.18~ld10,data=tAlla,family='binomial')
fitG1
BIC(fitG1)

nd$pGreen <- predict(fitG1,nd,type='response')
lines(pGreen~ld10,data=nd,lwd=2,col='green')

par(op)

#### NOW TRY MULTINOMIAL### SKIP THIS SECTION-###
require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(tAll[,c('ld10','PFstatus','Species.13','FireSev')],1,allNotNA)
table(rComp)
tAlla <- tAll[rComp,]

fit1 <- multinom(PFstatus ~ ld10,data=tAlla)
fit1
head(round(fitted(fit1),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit1)[,1],col='red')
points(tAlla$ld10,fitted(fit1)[,2],col='gray')
points(tAlla$ld10,fitted(fit1)[,3],col='green')


##########
head(tAll)
dim(tAll)
tAll$dupStatusCheck <- apply(tAll[,c('Dead.18','TopkillLive.18','gCrown.18')],1,sum,na.rm=T)
table(tAll$dupStatusCheck)
## same 2 as above

###end skip section ###

#### PROVISIONALLY ASSIGN TO THREE CLASSES
tAll$PFstatus <- (-1)
tAll$PFstatus[which(tAll$Live.18==0)] <- 0
tAll$PFstatus[which(tAll$TopkillLive.18==1)] <- 1
tAll$PFstatus[which(tAll$Topkill.18==0 & tAll$gCrown.18==1)] <- 2
table(tAll$PFstatus)

# now assign all remaing NAs to dead - TEMP STEP
tAll$PFstatus[which(tAll$PFstatus==(-1))] <- 0
table(tAll$PFstatus)

## NOW TRY MULTINOMIAL
# setup discrete FireSev
summary(tAll$FireSev)
tAll$dFS <- cut(tAll$FireSev,breaks = c(-200,35,130,298,1000))
table(tAll$dFS)

require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(tAll[,c('ld10','PFstatus','Species.13','FireSev')],1,allNotNA)
table(rComp)
tAlla <- tAll[rComp,]

fit1 <- multinom(PFstatus ~ ld10,data=tAlla)
fit1
head(round(fitted(fit1),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit1)[,1],col='red')
points(tAlla$ld10,fitted(fit1)[,2],col='gray')
points(tAlla$ld10,fitted(fit1)[,3],col='green')

#Size and Fire Sev 
fit2 <- multinom(PFstatus ~ ld10 + dFS,data=tAlla)
fit2
BIC(fit2)
head(round(fitted(fit2),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit2)[,1],col='red')
points(tAlla$ld10,fitted(fit2)[,2],col='grey')
points(tAlla$ld10,fitted(fit2)[,3],col='green')

# Add Species
fit3 <- multinom(PFstatus ~ ld10 + dFS + Species.13,data=tAlla)
fit3
BIC(fit3)
head(round(fitted(fit3),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit3)[,1],col='red')
points(tAlla$ld10,fitted(fit3)[,2],col='grey')
points(tAlla$ld10,fitted(fit3)[,3],col='green')

# Add Species * size interaction
fit3x <- multinom(PFstatus ~ ld10 + dFS + Species.13 + dFS:ld10,data=tAlla)
fit3x
BIC(fit3x)
head(round(fitted(fit3x),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit3x)[,1],col='red')
points(tAlla$ld10,fitted(fit3x)[,2],col='black')
points(tAlla$ld10,fitted(fit3x)[,3],col='green')

# COMPARE BIC
BIC(fit1)
BIC(fit2)
BIC(fit3)
BIC(fit3x)

# plot for selected values
rsel <- which(tAlla$Species.13=='PSEMEN' & as.numeric(tAlla$dFS)>0)
plot(tAlla$ld10[rsel],tAlla$PFstatus[rsel])
points(tAlla$ld10[rsel],fitted(fit3x)[rsel,1],col='red',pch=19)
points(tAlla$ld10[rsel],fitted(fit3x)[rsel,2],col='black',pch=19)
points(tAlla$ld10[rsel],fitted(fit3x)[rsel,3],col='green',pch=19)

## logit model and plot
d=tAll
xcn='ld10'
ycn='Live.18'
## NOT WORKING RIGHT NOW
logit2Plot <- function(d,xcn,ycn,np=101)
{
  dx <- d[,xcn]
  dy <- d[,ycn]
  plot(dy~dx)
  fit <- glm(dy~dx,family='binomial')
  nd <- data.frame(dx=seq(min(dx,na.rm=T),max(dx,na.rm=T),length.out=np))
  nd$yPred <- predict(fit,nd)
  lines(nd$dx,nd$yPred)
}
logit2Plot(tAll,'ld10','Live.18')



