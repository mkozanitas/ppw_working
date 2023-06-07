#analyzeData-tAll is to compare 2013, 2018, 2019, and 2020 for delayed mortality and other changes in state

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
newIndvs <- which(is.na(tAll$Plot.13))
length(newIndvs)

tAll$Plot.13[newIndvs] <- tAll$Plot.18[newIndvs]

# SEARCH AND REPLACE FROM HERE
tAll$Quad.13[newIndvs] <- tAll$Quad.18[newIndvs]
tAll$Type.13[newIndvs] <- tAll$Type.18[newIndvs]
tAll$Species.13[newIndvs] <- tAll$Species.18[newIndvs]

# This introduces inaccuracies because 2018 basal areas reflect 5 more years of growth
tAll$Basal.Area.13[newIndvs] <- tAll$Basal.Area.18[newIndvs]
tAll$dbh.13[newIndvs] <- tAll$dbh.18[newIndvs]
tAll$Dead.13[newIndvs] <- 0
tAll$Live.13[newIndvs] <- 1
tAll$gCrown.13[newIndvs] <- 1

# ignore these individuals for basal area growth
tAll$UseForBAGrowth <- T
tAll$UseForBAGrowth[newIndvs] <- F

## END CREATE PROXY 2103 Data for new trees encountered post-fire, including new plots


## quick look, out of curiosity, for basal area growth
#plot(tAll$Basal.Area.13[tAll$UseForBAGrowth],tAll$Basal.Area.18[tAll$UseForBAGrowth],log='xy')
#summary(tAll$Basal.Area.18/tAll$Basal.Area.13)
#abline(0,1,col='red')

## THIS CODE SECTION ANALYZES 2018 FATES ###

# ANALYZE BY SPECIES AND TYPE for 2018 post-fire fates
(use.species <- spNames$spName)

# create fst dataframe - FateSummaryTable
fst12 <- data.frame(Species=rep(use.species,each=2),Type=rep(c('SA','TR'),length(use.species)),N13=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
#fst12 <- data.frame(Species=rep(use.species,each=3),Type=rep(c('SA','TR','TS'),length(use.species)),N13=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
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
sum(fst12$nMissing) # 1330 duplicates - ignoring!!

# Summary across all species and types
tree.sum <- apply(fst12[fst12$Type=='TR',-c(1:2)],2,sum)
(tree.sum-tree.sum[6])/(tree.sum[1]-tree.sum[6])

sap.sum <- apply(fst12[fst12$Type=='SA',-c(1:2)],2,sum)
(sap.sum-sap.sum[6])/(sap.sum[1]-sap.sum[6])

#ts.sum <- apply(fst12[fst12$Type=='TS',-c(1:2)],2,sum)
#(ts.sum-sap.sum[6])/(ts.sum[1]-ts.sum[6])

fate.sum <- apply(fst12[,-c(1:2)],2,sum)
(fate.sum-fate.sum[6])/(fate.sum[1]-fate.sum[6])

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

# plots experiencing each fire severity level, and how many N13 individuals in each
table(tAll$Plot.13[which(tAll$fsLevel==0)])
table(tAll$Plot.13[which(tAll$fsLevel==1)])
table(tAll$Plot.13[which(tAll$fsLevel==2)])
table(tAll$Plot.13[which(tAll$fsLevel==3)])

# Use dbh for largest stem
nodbh <- which(is.na(tAll$dbh.13))
length(nodbh)
head(tAll[nodbh,])
summary(tAll$dbh.13)

## USE LOG DBH for analysis
hist(tAll$dbh.13)
tAll$ldbh <- log10(tAll$dbh.13)
hist(tAll$ldbh)
tAll$ldbh2 <- tAll$ldbh^2

table(tAll$Live.18, tAll$fate.18)

plot(tAll$ldbh,tAll$Live.18)
which(tAll$Live.18==0&tAll$fate.18=="LN")
# all clear

#### SURVIVAL ANALYSIS
(N <- length(which(!is.na(tAll$ldbh) & !is.na(tAll$Live.18))))
fit <- glm(Live.18~ldbh,data=tAll,family='binomial')
summary(fit)

summary(tAll$ldbh)
nd <- with(tAll,data.frame(ldbh=seq(min(tAll$ldbh,na.rm=T),max(tAll$ldbh,na.rm=T),length.out=1000)))
head(nd)
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)

range(tAll$ldbh,na.rm=T)
plot(tAll$ldbh,tAll$Live.18)
lines(nd$ldbh,nd$pSurvAll,lwd=4)
#abline(h=0.5,lty=2) 
#h draws horizontal line at .5, lty -dashed or solid, lwd is line width

#what is critical basal area to achieve 50% survival?
(ld50 <- nd$ldbh[which(nd$pSurvAll>=0.5)[1]])
abline(v=ld50,lty=2)

(spRes <- data.frame(species='All',N=N,ld50=ld50,slp=fit$coefficients[2]))
spRes

## now run by species for species with lots of data
spN <- table(tAll$Species.13)
spN[order(spN)]

(spA <- names(spN)[which(spN>=25)])
spA <- spA[-which(spA=='QUEBER')] # drop QUEBER - none died!

## Abundant species only
#AbSp <- c('AMOCAL','ARBMEN','ARCMAN','FRACAL','HETARB','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')

# dropping QUEBER, FRACAL, AMOCAL for now
AbSp <- c('ARBMEN','ARCMAN','HETARB','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')
length(AbSp)

# use abundant species identified above for ESA
(spA <- AbSp)

i <- 3
for (i in 1:length(spA))
{
  (species <- spA[i])
  tmp <- tAll[which(tAll$Species.13==species),]
  N <- length(which(!is.na(tmp$Basal.Area.13) & !is.na(tmp$Live.18)))
  fit <- glm(Live.18~ldbh,data=tmp,family='binomial')
  pval <- predict(fit,newdata=nd,type='response')
  nd <- data.frame(nd,pval)
  names(nd)[length(names(nd))] <- paste('pSurv_',species,sep='')
  lines(nd$ldbh,nd[,ncol(nd)])
  
  ld50 <- NA
  if (fit$coefficients[2]>0) 
    if (nd[1,ncol(nd)]<0.5) 
      ld50 <- nd$ldbh[which(nd[,ncol(nd)]>=0.5)[1]] 
  spRes <- rbind(spRes,c(species,N,ld50,fit$coefficients[2]))
}
spRes$ld50 <- as.numeric(spRes$ld50)
spRes$d50 <- round(10^spRes$ld50,3)
spRes$ld50 <- round(as.numeric(spRes$ld50),3)
spRes$slp <- round(as.numeric(spRes$slp),3)
spRes

plotSP <- function(tAll,species=NULL,xlims=range(tAll$ldbh,na.rm=T))
{
  tmp <- tAll[which(tAll$Species.13==species),]
  mindbh <- min(tmp$ldbh,na.rm=T)
  maxdbh <- max(tmp$ldbh,na.rm=T)
  plot(tmp$ldbh,tmp$Live.18,main=species,xlim=xlims)
  fit <- glm(Live.18~ldbh,data=tmp,family='binomial')
  nd <- with(tAll,data.frame(ldbh=seq(mindbh,maxdbh,length.out=1001)))
  pval <- predict(fit,newdata=nd,type='response')
  lines(nd$ldbh,pval)  
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

########### Starting with Models########

## Full model with species (can subset out data with mod/high fslevel(2:3) for visualization only- FS not in model-change back to c(0:3 to expand to all fslevels) to generate curves
tAlls <- tAll[which(tAll$Species.13 %in% spA & tAll$fsLevel %in% c(2:3)),]

# assign dependent variable to rVar
names(tAlls)
selVar <- 'Live.18'
tAlls$rVar <- tAlls[,selVar]

(N <- length(which(!is.na(tAlls$ldbh) & !is.na(tAlls$rVar))))
table(tAlls$Species.13)
head(tAlls)
fit <- glm(rVar~ldbh + ldbh2 + Species.13,data=tAlls,family='binomial')
summary(fit)

summary(tAlls$ldbh)
#choose sizes here at which predicted values will be calculated
predSizes <- c(10,30,80)
nd <- with(tAlls,data.frame(ldbh=rep(log10(predSizes),length(spA)),Species.13=rep(spA,each=length(predSizes))))
nd$ldbh2 <- nd$ldbh^2
nd
nd$pValue <- round(predict(fit,newdata=nd,type='response'),4)
#- then use Sp.names instead of Species.13- below to reorder for barplot
nd$Sp.Names <- c("aARCMAN","bPSEMEN", "cQUEDOU", "dQUEKEL", "eARBMEN", "fQUEGAR", "gUMBCAL", "hHETARB", "iQUEAGR")


# with shrubs
selSize <- 3
barplot(pValue~Species.13,data=nd[which(nd$ldbh==log10(predSizes[selSize])),],ylim=c(0,1),main=paste(selVar,'predicted value at ',predSizes[selSize],' cm dbh'))

# remove shrubs for larger sizes cm plot
# op=par will stack 3 figures (1,3) for horizontal & (3,1) for vertical
op=par(mfrow=c(1,3)) 
selSize <- 1
barplot(pValue~Species.13,data=nd[which(nd$ldbh==log10(predSizes[selSize])  & !nd$Species.13 %in% c('HETARB','ARCMAN')),],ylim=c(0,1),main=paste(selVar,'predicted value at ',predSizes[selSize],' cm dbh'))

selSize <- 2
barplot(pValue~Species.13,data=nd[which(nd$ldbh==log10(predSizes[selSize])  & !nd$Species.13 %in% c('HETARB','ARCMAN')),],ylim=c(0,1),main=paste(selVar,'predicted value at ',predSizes[selSize],' cm dbh'))

selSize <- 3
barplot(pValue~Species.13,data=nd[which(nd$ldbh==log10(predSizes[selSize])  & !nd$Species.13 %in% c('HETARB','ARCMAN')),],ylim=c(0,1),main=paste(selVar,'predicted value at ',predSizes[selSize],' cm dbh'))


#run to reset to single figure 
par(op)


### MODELS WITH FIRE SEVERITY - use 4 levels as factors
# Change response variable here and then run model. Swap gCrown.18 or Live.18 to analyze crown survival
tAlls <- tAll[which(tAll$Species.13 %in% spA & tAll$fsLevel %in% c(0:3)),]
#did a fix here to set variable like we did above using selVar..
selVar <- 'gCrown.18'
tAlls$rVar <- tAlls[,selVar]
fit <- glm(rVar~ldbh + as.factor(fsLevel),data=tAlls,family='binomial')
summary(fit)

nvals <- 101
ldbh.vals <- seq(-0.5,2,length.out=nvals)
#FireSev.vals <- seq(min(fs[,fsmet],na.rm=T),max(fs[,fsmet],na.rm=T),length.out=nvals)
FireLevel.vals <- sort(unique(tAlls$fsLevel))
nd <- with(tAll,data.frame(ldbh=rep(ldbh.vals,length(FireLevel.vals)),fsLevel=rep(FireLevel.vals,each=nvals)))
head(nd)

nd$rVarPred <- predict(fit,newdata=nd,type='response')
head(nd)
tail(nd)

# Plot response variable as a function of size, with isoclines as a function fire severity. Lines are increasing - survival is higher for larger trees, but lower at higher fire severity
plot(tAlls$ldbh,tAlls$rVar,type="n", xlab="log10DBH", ylab = selVar)
i=1
for (i in 1:length(FireLevel.vals)) {
  ndt <- nd[which(nd$fsLevel==FireLevel.vals[i]),]
  lines(ndt$ldbh,ndt$rVarPred,col="red")
  text(ldbh.vals[50],ndt$rVarPred[which(ndt$ldbh==ldbh.vals[50])],FireLevel.vals[i])
}

#### THIS WON'T WORK NOW AS MODEL ABOVE WAS CHANGED TO USE FIRE SEVERITY LEVELS
# Plot survival as a function of fire severity, with isoclines as a function ldbh. Lines are declining - survival is lower at higher fire severity, but larger trees have higher values
# plot(tAlls$FireSev,tAlls$Live.18)
# i=1
# for (i in 1:nvals) {
#   ndt <- nd[which(nd$ldbh==ldbh.vals[i]),]
#   lines(ndt$FireSev,ndt$pSurvAll)
#   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
# }
### END COMMENT OUT

## RECODED WITH FIRE LEVELS
# model with size, fire, species, predicting survival
{
  fit2 <- glm(Live.18~ldbh + as.factor(fsLevel) + Species.13,data=tAlls,family='binomial')
  fit1 <- glm(Live.18~ldbh + as.factor(fsLevel) + Species.13+ as.factor(fsLevel):Species.13,data=tAlls,family='binomial')
  BIC(fit1)
  BIC(fit2)
  summary(fit1)
  summary(fit2)
  
  nvals <- 11
  ldbh.vals <- seq(-0.5,2,length.out=nvals)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  plotSpecies <- function(spname,nd.tmp=nd2) {
    tmp <- tAlls[which(tAlls$Species.13==spname),]
    #op=par(mfrow=c(1,2))
    
    # plot(tmp$FireSev,tmp$Live.18,main=spname,xlim=range(tAlls$FireSev,na.rm=T))
    # i=1
    # for (i in 1:nvals) {
    #   ndt <- nd2[which(nd2$Species.13==spname & nd2$ldbh==ldbh.vals[i]),]
    #   lines(ndt$FireSev,ndt$pSurvAll)
    #   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
    # }
    
    plot(tmp$ldbh,tmp$Live.18,xlim=range(tAlls$ldbh,na.rm=T))
    i=2
    for (i in 1:length(unique(nd.tmp$fsLevel))) {
      ndt <- nd.tmp[which(nd2$Species.13==spname & nd.tmp$fsLevel==fsLevels.vals[i]),]
      lines(ndt$ldbh,ndt$pSurvAll)
      text(ldbh.vals[5],ndt$pSurvAll[which(ndt$ldbh==ldbh.vals[5])],fsLevels.vals[i])
    }
    #par(op)
  }
  #This will plot isoclines of survival at each fire severity for any indv species 
  plotSpecies('PSEMEN')
  plotSpecies('ARBMEN')
  plotSpecies('QUEAGR')
  plotSpecies('UMBCAL')
  
  #using data from model w/ fire sev included- not subsetted out for visualization as above
  # obtain predicted main effects of fire severity for each species at a common size
  nvals <- 1
  ldbh.vals <- log10(2)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
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
  fit2 <- glm(gCrown.18~ldbh + as.factor(fsLevel) + Species.13,data=tAlls,family='binomial')
  fit1 <- glm(gCrown.18~ldbh + as.factor(fsLevel) + Species.13+ as.factor(fsLevel):Species.13,data=tAlls,family='binomial')
  BIC(fit1)
  BIC(fit2)
  summary(fit1)
  summary(fit2)
  
  nvals <- 11
  ldbh.vals <- seq(-0.5,2,length.out=nvals)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  spname='QUEAGR'
  plotSpecies <- function(spname,nd.tmp=nd2) {
    tmp <- tAlls[which(tAlls$Species.13==spname),]
    #op=par(mfrow=c(1,2))
    
    # plot(tmp$FireSev,tmp$Live.18,main=spname,xlim=range(tAlls$FireSev,na.rm=T))
    # i=1
    # for (i in 1:nvals) {
    #   ndt <- nd2[which(nd2$Species.13==spname & nd2$ldbh==ldbh.vals[i]),]
    #   lines(ndt$FireSev,ndt$pSurvAll)
    #   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
    # }
    
    plot(tmp$ldbh,tmp$gCrown.18,xlim=range(tAlls$ldbh,na.rm=T))
    i=2
    for (i in 1:length(unique(nd.tmp$fsLevel))) {
      ndt <- nd.tmp[which(nd2$Species.13==spname & nd.tmp$fsLevel==fsLevels.vals[i]),]
      lines(ndt$ldbh,ndt$pSurvAll)
      text(ldbh.vals[5],ndt$pSurvAll[which(ndt$ldbh==ldbh.vals[5])],fsLevels.vals[i])
    }
    #par(op)
  }
  plotSpecies('QUEGAR')
  
  # obtain predicted main effects of fire severity for each species at a common size (2CM)
  nvals <- 1
  ldbh.vals <- log10(2)
  fsLevels.vals <- c(0:3)
  nd <- with(tAll,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.13=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
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


(PFstatusLevels <- 0:2)
(PFsPlotVals <- c(0.95,1,1.05))
(PFsPlotCols <- c('black','red','green'))

tAll$PFsPlotVals <- PFsPlotVals[match(tAll$PFstatus,PFstatusLevels)]
tAll$PFsPlotCols <- PFsPlotCols[match(tAll$PFstatus,PFstatusLevels)]

# range of ldbh for abundant species
spArows <- which(tAll$Species.13 %in% spA)
(tAllldbh.range <- c(min(tAll$ldbh[spArows],na.rm=T),max(tAll$ldbh[spArows],na.rm=T)))

# head(tAll$fsLevel)
# head(match(tAll$fsLevel,PFstatusLevels))
# head(tAll$PFsPlotVals)
## ANALYSIS FOR ONE SPECIES
spA

# pick one of these!
selSpecies <- spA # use spA for all abundant species, rather than one species
FireLevels <- c('Mod+High')
#FireLevels <- c('ANY LEVEL') #then change range in line below c(1:3)
tAllsp <- tAll[which(tAll$Species.13 %in% c(selSpecies) & tAll$fsLevel %in% c(2:3)),] #individual species?
{
  #tAllsp <- tAll[which(tAll$Species.13 %in% spA),] #abundant species?
  #tAllsp <- tAll[which(tAll$Species.13 %in% spA & tAll$fsLevel>1),] #abundant sp with fs level of 1 or more?
  
  dim(tAllsp)
  tAllsp <- tAllsp[which(!is.na(tAllsp$ldbh)),]
  dim(tAllsp)
  
  # code to run a binomial model and plot response curve with data
  # MORTALITY
  nd <- with(tAllsp,data.frame(ldbh=seq(min(tAllsp$ldbh,na.rm=T),max(tAllsp$ldbh,na.rm=T),length.out=101)))
  nd$ldbh2 <- nd$ldbh^2
  
  #plot(tAllsp$ldbh,tAllsp$Dead.18,xlim=tAllldbh.range)
  fit1 <- glm(Dead.18~ldbh,data=tAllsp,family='binomial')
  fit2 <- glm(Dead.18~ldbh+ldbh2,data=tAllsp,family='binomial')
  BIC(fit1)
  BIC(fit2)
  nd$pMortality <- predict(fit1,newdata=nd,type='response')
  #lines(nd$ldbh,nd$pMortality)
  
  #TOPKILL WITH RESPROUT
  #plot(tAllsp$ldbh,tAllsp$TB,xlim=tAllldbh.range)
  fit <- glm(DR.18~ldbh+ldbh2 ,data=tAllsp,family='binomial')
  summary(fit)
  nd$pDR <- predict(fit,newdata=nd,type='response')
  #lines(nd$ldbh,nd$pTB)
  
  #GREEN CROWN
  #plot(tAllsp$ldbh,tAllsp$gCrown.18,xlim=tAllldbh.range)
  fit1 <- glm(gCrown.18~ldbh,data=tAllsp,family='binomial')
  fit2 <- glm(gCrown.18~ldbh+ldbh2,data=tAllsp,family='binomial')
  summary(fit1)
  nd$pGCrown <- predict(fit1,newdata=nd,type='response')
  #lines(nd$ldbh,nd$pGCrown)
  
  #plot all three (change main from selSpecies to ".." to alter main title)
  plot(tAllsp$ldbh,tAllsp$PFsPlotVals,col=tAllsp$PFsPlotCols,pch=19,ylim=c(-0.05,1.05),xlim=tAllldbh.range,main=paste("Allsp",FireLevels))
  points(tAllsp$ldbh,rep(-0.05,length(tAllsp$ldbh)))
  lines(nd$ldbh,nd$pGCrown,col='green')
  lines(nd$ldbh,nd$pDR,col='red')
  lines(nd$ldbh,nd$pMortality)
  
  # does the sum of the three binomials for these three exclusive fates sum to 1?
  nd$pTOT <- apply(nd[,c('pGCrown','pDR','pMortality')],1,sum)
  summary(nd$pTOT)
}

#reset plotting window with one panel
par(mfrow=c(1,1))
# NOW FIT MULTINOMIAL
require(nnet)

# MULTINOMIAL - QUADRATIC CAN BE ADDED HERE '+ldbh2' - changes results some
fit1 <- multinom(as.factor(PFstatus.18) ~ ldbh +ldbh2, data=tAllsp)
fit1
head(round(fitted(fit1),2))
dim(fitted(fit1))

plot(tAllsp$ldbh,tAllsp$Live.18)
points(tAllsp$ldbh,fitted(fit1)[,1],col='black')
points(tAllsp$ldbh,fitted(fit1)[,2],col='red')
points(tAllsp$ldbh,fitted(fit1)[,3],col='green')
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
#fit <- glmer(Live.18~ldbh + FireSev + Species.13 + (1 |Plot.13),data=tAlls,family='binomial')

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
tAlla$ldbh2 <- tAlla$ldbh^2

op=par(mfrow=c(1,1))
plot(Dead.18~ldbh,data=tAlla)
fitD <- glm(Dead.18~ldbh+ldbh2,data=tAlla,family='binomial')
fitD
nd <- data.frame(ldbh=seq(min(tAlla$ldbh,na.rm=T),max(tAlla$ldbh,na.rm=T),length.out=101))
nd$ldbh2 <- nd$ldbh^2
head(nd)
nd$pDead <- predict(fitD,nd,type='response')
lines(pDead~ldbh,data=nd,lwd=2,col='black')

#plot(Topkill.18~ldbh,data=tAlla)
fitT <- glm(TopkillLive.18~ldbh+ldbh2,data=tAlla,family='binomial')
fitT
nd$pTopKill <- predict(fitT,nd,type='response')
lines(pTopKill~ldbh,data=nd,lwd=2,col='red')

#plot(gCrown.18~ldbh,data=tAlla)
fitG2 <- glm(gCrown.18~ldbh+ldbh2,data=tAlla,family='binomial')
fitG2
BIC(fitG2)
fitG1 <- glm(gCrown.18~ldbh,data=tAlla,family='binomial')
fitG1
BIC(fitG1)

nd$pGreen <- predict(fitG1,nd,type='response')
lines(pGreen~ldbh,data=nd,lwd=2,col='green')

par(op)

#### NOW TRY MULTINOMIAL### SKIP THIS SECTION-###
require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(tAll[,c('ldbh','PFstatus','Species.13','FireSev')],1,allNotNA)
table(rComp)
tAlla <- tAll[rComp,]

fit1 <- multinom(PFstatus ~ ldbh,data=tAlla)
fit1
head(round(fitted(fit1),2))

plot(tAlla$ldbh,tAlla$Live.18)
points(tAlla$ldbh,fitted(fit1)[,1],col='red')
points(tAlla$ldbh,fitted(fit1)[,2],col='gray')
points(tAlla$ldbh,fitted(fit1)[,3],col='green')


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
rComp <- apply(tAll[,c('ldbh','PFstatus','Species.13','FireSev')],1,allNotNA)
table(rComp)
tAlla <- tAll[rComp,]

fit1 <- multinom(PFstatus ~ ldbh,data=tAlla)
fit1
head(round(fitted(fit1),2))

plot(tAlla$ldbh,tAlla$Live.18)
points(tAlla$ldbh,fitted(fit1)[,1],col='red')
points(tAlla$ldbh,fitted(fit1)[,2],col='gray')
points(tAlla$ldbh,fitted(fit1)[,3],col='green')

#Size and Fire Sev 
fit2 <- multinom(PFstatus ~ ldbh + dFS,data=tAlla)
fit2
BIC(fit2)
head(round(fitted(fit2),2))

plot(tAlla$ldbh,tAlla$Live.18)
points(tAlla$ldbh,fitted(fit2)[,1],col='red')
points(tAlla$ldbh,fitted(fit2)[,2],col='grey')
points(tAlla$ldbh,fitted(fit2)[,3],col='green')

# Add Species
fit3 <- multinom(PFstatus ~ ldbh + dFS + Species.13,data=tAlla)
fit3
BIC(fit3)
head(round(fitted(fit3),2))

plot(tAlla$ldbh,tAlla$Live.18)
points(tAlla$ldbh,fitted(fit3)[,1],col='red')
points(tAlla$ldbh,fitted(fit3)[,2],col='grey')
points(tAlla$ldbh,fitted(fit3)[,3],col='green')

# Add Species * size interaction
fit3x <- multinom(PFstatus ~ ldbh + dFS + Species.13 + dFS:ldbh,data=tAlla)
fit3x
BIC(fit3x)
head(round(fitted(fit3x),2))

plot(tAlla$ldbh,tAlla$Live.18)
points(tAlla$ldbh,fitted(fit3x)[,1],col='red')
points(tAlla$ldbh,fitted(fit3x)[,2],col='black')
points(tAlla$ldbh,fitted(fit3x)[,3],col='green')

# COMPARE BIC
BIC(fit1)
BIC(fit2)
BIC(fit3)
BIC(fit3x)

# plot for selected values
rsel <- which(tAlla$Species.13=='PSEMEN' & as.numeric(tAlla$dFS)>0)
plot(tAlla$ldbh[rsel],tAlla$PFstatus[rsel])
points(tAlla$ldbh[rsel],fitted(fit3x)[rsel,1],col='red',pch=19)
points(tAlla$ldbh[rsel],fitted(fit3x)[rsel,2],col='black',pch=19)
points(tAlla$ldbh[rsel],fitted(fit3x)[rsel,3],col='green',pch=19)

## logit model and plot
d=tAll
xcn='ldbh'
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
logit2Plot(tAll,'ldbh','Live.18')



