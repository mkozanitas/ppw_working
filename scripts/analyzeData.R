# RUN prepareData and examineData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R

rm(list=ls())
library(lme4)

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

all.id <- readRDS('data/allid-nodups.Rdata')
str(all.id)

spNames <- read.csv('data/all-spp-names.csv')
head(spNames)
#allIndv <- readRDS('data/allIndv.Rdata')
allIndv <- read.csv('data/allIndv.csv')
head(allIndv)

#### First round of analysis of post-fire fates, 2013-2018
## get percent survival by species and type

# select which individuals to eliminate - this can be done separately for each analysis
# take any individuals where plot and species match for the first two years
nrow(allIndv)
plot.ok <- allIndv$Num[which(allIndv$P13==allIndv$P18)] 
spec.ok <- allIndv$Num[which(allIndv$S13==allIndv$S18)] 
ps.ok <- intersect(plot.ok,spec.ok)
length(ps.ok)

# include trees from new plots
ps.ok <- c(ps.ok,allIndv$Num[which(allIndv$P18 %in% c('PPW1851','PPW1852','PPW1853','PPW1854'))])
length(ps.ok)

head(ps.ok)
tail(ps.ok)

# were any individuals missed in 2018 included in list above? (NA18=1 if Num is present in 13 and 19, and missing in 18)
na18 <- all.id[[1]]$Num[which(all.id[[1]]$NA18==1)]
length(na18)
na18[which(na18 %in% ps.ok)]
# none of them captured in ps.ok, so we can ignore them for this analysis

t1 <- all.id[[1]][which(all.id[[1]]$Num %in% ps.ok),]
nrow(t1)
t2 <- all.id[[2]][which(all.id[[2]]$Num %in% ps.ok),]
nrow(t2)

# And merge!
t12 <- merge(t1,t2,by = 'Num',all = T)
dim(t12)
head(t12)
tail(t12)

# create 'fake 2013 data'
# rownums for new 2018 individuals from new plots
newIndvs <- which(is.na(t12$Plot.x))
length(newIndvs)

t12[4042:4043,]
t12$Plot.x[newIndvs] <- t12$Plot.y[newIndvs]
t12$Quad.x[newIndvs] <- t12$Quad.y[newIndvs]
t12$Type.x[newIndvs] <- t12$Type.y[newIndvs]
t12$Species.x[newIndvs] <- t12$Species.y[newIndvs]

# This introduces error because 2018 basal areas reflect 5 more years of growth
t12$Basal.Area.x[newIndvs] <- t12$Basal.Area.y[newIndvs]
t12$Dead.x[newIndvs] <- 0
t12$Live.x[newIndvs] <- 1
t12$gCrown.x[newIndvs] <- 1

#plot(t12$Basal.Area.x,t12$Basal.Area.y)
#summary(t12$Basal.Area.y/t12$Basal.Area.x)
#abline(0,1)

# ANALYZE BY SPECIES AND TYPE
(use.species <- spNames$x)
fst <- data.frame(Species=rep(use.species,each=2),Type=rep(c('SA','TR'),length(use.species)),N13=NA,N18.dead=NA,N18.TKB=NA,N18.BC=NA,N18.C=NA,nMissing=NA)
head(fst)                  
tail(fst)

i=24
for (i in 1:nrow(fst))
{
  sp <- fst$Species[i]
  ty <- fst$Type[i]
  init <- t12$Num[which(t12$Species.x==sp & t12$Type.x==ty & t12$Live.x==1)]
  fst$N13[i] <- length(init)
  
  fin1 <- intersect(init,t12$Num[which(t12$Live.y==0)])
  # fin1 is the number of original ty surviving, whether or not they transitioned from SA->TR (or the other way!)
  
  fst$N18.dead[i] <- length(fin1)
  
  #fin2 is topkill and basal sprouting
  fin2 <- intersect(init,t12$Num[which(t12$bSprout.y==1 & t12$Topkill.y==1)])
  fst$N18.TKB[i] <- length(fin2) 
  
  # fin3 is basal sprouting and green crown
  fin3 <- intersect(init,t2$Num[which(t2$bSprout==1 & t2$gCrown==1)])
  fst$N18.BC[i] <- length(fin3) 
  
  # fin4 is green crown only
  fin4 <- intersect(init,t2$Num[which(t2$bSprout==0 & t2$gCrown==1)])
  fst$N18.C[i] <- length(fin4) 
  
  missed <- init[which (!init %in% c(fin1,fin2,fin3,fin4))]
  fst$nMissing[i] <- length(missed)
}
head(fst)
tail(fst)
sum(fst$nMissing)

fst$percSurv <- 1 - fst$N18.dead/fst$N13
fst$percSurv[fst$N13==0] <- NA

## What percent died in Tubbs, overall and by Type?
SArows <- which(fst$Type=='SA')
TRrows <- which(fst$Type=='TR')

sum(fst$N13)
sum(fst$N18.dead)
sum(fst$N18.dead)/sum(fst$N13)
sum(fst$N18.dead[SArows])/sum(fst$N13[SArows])
sum(fst$N18.dead[TRrows])/sum(fst$N13[TRrows])

## Abundant species only
#AbSp <- c('AMOCAL','ARBMEN','ARCMAN','FRACAL','HETARB','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')

# dropping QUEBER, FRACAL, AMOCAL for now
AbSp <- c('ARBMEN','ARCMAN','HETARB','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')
length(AbSp)
fsta <- fst[which(fst$Species %in% AbSp),]

#examine data
fsta[fsta$Type=='TR',][order(fsta$percSurv[fsta$Type=='TR']),]
fsta[fsta$Type=='SA',][order(fsta$percSurv[fsta$Type=='SA']),]

# plot sapling vs. tree survival
plot(fsta[fsta$Type=='TR','percSurv'],fsta[fsta$Type=='SA','percSurv'],xlim=c(0,1),ylim=c(0,1),type='n',xlab='Survival, trees',ylab='Survival, saplings')
text(fsta[fsta$Type=='TR','percSurv'],fsta[fsta$Type=='SA','percSurv'],labels = fsta[fsta$Type=='TR','Species'])
abline(0,1)

#### EXPAND ANALYSIS to ALL POST-FIRE FATES
dim(t12)
head(t12)
names(t12)

## choose fire severity metric
fsmet <- 'Tubbs.MTBS.RDNBR.30'
names(fs)
summary(fs[,fsmet])
hist(fs[,fsmet])
sort(fs[,fsmet])

f2t <- match(t12$Plot.x,fs$Plot)
head(f2t)
tail(f2t)
t12$FireSev <- fs[f2t,fsmet]
dim(t12)
tail(t12)
table(t12$Plot.x)

# RDNBR fire severity levels
# Unburned
# Low <170, but burned
# Medium 170-700
# High: >700
# unburned: 1308, 1309, 1311, 1312, 1327, 1344, 1347

summary(t12$FireSev)
hist(t12$FireSev,breaks=c(-100,0,50,100,150,200,300,400,500,600,700,800,1000))
fs.breaks <- c(-100,135,430,1000) # Based on Parks et al. 2014, but not using intermediate split at 304
t12$fsCat <- cut(t12$FireSev,fs.breaks)
table(t12$fsCat)
t12$fsLevel <- as.numeric(t12$fsCat)
table(t12$fsLevel)

# manually code unburned as 0
unburned.plots <- c('PPW1308','PPW1309','PPW1311','PPW1312','PPW1327','PPW1344','PPW1347')
t12$fsLevel[which(t12$Plot.x %in% unburned.plots)] <- 0
table(t12$fsLevel)

table(t12$Plot.x[which(t12$fsLevel==3)])

# CALC DBH - easier to connect with our field data and knowledge
t12$dbh <- ba2d(t12$Basal.Area.x)
summary(t12$dbh)

## USE LOG DBH for analysis
hist(t12$dbh)
t12$ldbh <- log10(t12$dbh)
hist(t12$ldbh)
t12$ldbh2 <- t12$ldbh^2

#### SURVIVAL ANALYSIS
plot(t12$ldbh,t12$Live.y)

# **** mixing SA and TR - need to work on diameter conversion to get this right ***
(N <- length(which(!is.na(t12$ldbh) & !is.na(t12$Live.y))))
fit <- glm(Live.y~ldbh,data=t12,family='binomial')
summary(fit)

summary(t12$ldbh)
nd <- with(t12,data.frame(ldbh=seq(min(t12$ldbh,na.rm=T),max(t12$ldbh,na.rm=T),length.out=1000)))
head(nd)
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)

plot(t12$ldbh,t12$Live.y)
lines(nd$ldbh,nd$pSurvAll,lwd=4)
abline(h=0.5,lty=2)

#what is critical basal area to achieve 50% survival?
(ld50 <- nd$ldbh[which(nd$pSurvAll>=0.5)[1]])
abline(v=ld50,lty=2)

(spRes <- data.frame(species='All',N=N,ld50=ld50,slp=fit$coefficients[2]))
spRes

## now run by species for species with lots of data
spN <- table(t12$Species.x)
spN[order(spN)]

#(spA <- names(spN)[which(spN>=25)])
#spA <- spA[-which(spA=='QUEBER')] # drop QUEBER - none died!

# use abundant species identified above for ESA
(spA <- AbSp)

i <- 3
for (i in 1:length(spA))
{
  (species <- spA[i])
  tmp <- t12[which(t12$Species.x==species),]
  N <- length(which(!is.na(tmp$Basal.Area.x) & !is.na(tmp$Live.y)))
  fit <- glm(Live.y~ldbh,data=tmp,family='binomial')
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

plotSP <- function(t12,species=NULL)
{
  tmp <- t12[which(t12$Species.x==species),]
  plot(tmp$ldbh,tmp$Live.y)
  fit <- glm(Live.y~ldbh,data=tmp,family='binomial')
  nd <- with(t12,data.frame(ldbh=seq(-0.5,2,length.out=1001)))
  pval <- predict(fit,newdata=nd,type='response')
  lines(nd$ldbh,pval)  
}
plotSP(t12,'ARCMAN')
plotSP(t12,'QUEAGR')
plotSP(t12,'QUEDOU')
plotSP(t12,'ARBMEN')
plotSP(t12,'HETARB')
plotSP(t12,'PSEMEN')
plotSP(t12,'QUEGAR')
plotSP(t12,'QUEKEL')
plotSP(t12,'UMBCAL')

## Full model with species
t12s <- t12[which(t12$Species.x %in% spA),]
(N <- length(which(!is.na(t12s$ldbh) & !is.na(t12s$Live.y))))
head(t12s)

fit <- glm(Live.y~ldbh * Species.x,data=t12s,family='binomial')
summary(fit)

#nd <- with(t12s,data.frame(ldbh=seq(min(t12s$ldbh,na.rm=T),max(t12s$ldbh,na.rm=T),length.out=101)))
summary(t12s$ldbh)

nd <- with(t12s,data.frame(ldbh=log10(2),Species.x=sort(unique(Species.x))))
nd
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
##nd$Sp.Names <- c("A.ARCMAN")- then use Sp.names
barplot(nd$pSurvAll~nd$Species.x,ylim=c(0,0.8))

### MODELS WITH FIRE SEVERITY - use 4 levels as factors
# Change response variable here and then run model. Swap gCrown.y or Live.y to analyze crown survival
#t12s$rVar <- t12s$Live.y 
t12s$rVar <- t12s$gCrown.y 
fit <- glm(rVar~ldbh + as.factor(fsLevel),data=t12s,family='binomial')
summary(fit)

nvals <- 101
ldbh.vals <- seq(-0.5,2,length.out=nvals)
#FireSev.vals <- seq(min(fs[,fsmet],na.rm=T),max(fs[,fsmet],na.rm=T),length.out=nvals)
FireLevel.vals <- sort(unique(t12s$fsLevel))
nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,length(FireLevel.vals)),fsLevel=rep(FireLevel.vals,each=nvals)))
head(nd)

nd$rVarPred <- predict(fit,newdata=nd,type='response')
head(nd)
tail(nd)

# Plot response variable as a function of size, with isoclines as a function fire severity. Lines are increasing - survival is higher for larger trees, but lower at higher fire severity
plot(t12s$ldbh,t12s$rVar,type="n", xlab="log10DBH", ylab = "CrownSurvival")
i=1
for (i in 1:length(FireLevel.vals)) {
  ndt <- nd[which(nd$fsLevel==FireLevel.vals[i]),]
  lines(ndt$ldbh,ndt$rVarPred,col="red")
  text(ldbh.vals[50],ndt$rVarPred[which(ndt$ldbh==ldbh.vals[50])],FireLevel.vals[i])
}

#### THIS WON'T WORK NOW AS MODEL ABOVE WAS CHANGED TO USE FIRE SEVERITY LEVELS
# Plot survival as a function of fire severity, with isoclines as a function ldbh. Lines are declining - survival is lower at higher fire severity, but larger trees have higher values
# plot(t12s$FireSev,t12s$Live.y)
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
fit2 <- glm(Live.y~ldbh + as.factor(fsLevel) + Species.x,data=t12s,family='binomial')
fit1 <- glm(Live.y~ldbh + as.factor(fsLevel) + Species.x+ as.factor(fsLevel):Species.x,data=t12s,family='binomial')
BIC(fit1)
BIC(fit2)
summary(fit1)
summary(fit2)

nvals <- 11
ldbh.vals <- seq(-0.5,2,length.out=nvals)
fsLevels.vals <- c(0:3)
nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
dim(nd)
head(nd)

nd2 <- data.frame(Species.x=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))

nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
head(nd2)

spname='QUEAGR'
plotSpecies <- function(spname,nd.tmp=nd2) {
  tmp <- t12s[which(t12s$Species.x==spname),]
  #op=par(mfrow=c(1,2))
  
  # plot(tmp$FireSev,tmp$Live.y,main=spname,xlim=range(t12s$FireSev,na.rm=T))
  # i=1
  # for (i in 1:nvals) {
  #   ndt <- nd2[which(nd2$Species.x==spname & nd2$ldbh==ldbh.vals[i]),]
  #   lines(ndt$FireSev,ndt$pSurvAll)
  #   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
  # }
  
  plot(tmp$ldbh,tmp$Live.y,xlim=range(t12s$ldbh,na.rm=T))
  i=2
  for (i in 1:length(unique(nd.tmp$fsLevel))) {
    ndt <- nd.tmp[which(nd2$Species.x==spname & nd.tmp$fsLevel==fsLevels.vals[i]),]
    lines(ndt$ldbh,ndt$pSurvAll)
    text(ldbh.vals[5],ndt$pSurvAll[which(ndt$ldbh==ldbh.vals[5])],fsLevels.vals[i])
  }
  #par(op)
}
plotSpecies('QUEGAR')

# obtain predicted main effects of fire severity for each species at a common size
nvals <- 1
ldbh.vals <- log10(2)
fsLevels.vals <- c(0:3)
nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
dim(nd)
head(nd)

nd2 <- data.frame(Species.x=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))

nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
head(nd2)

barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==0),],ylim=c(0,1))
barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==1),],ylim=c(0,1))
barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==2),],ylim=c(0,1))
barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==3),],ylim=c(0,1))
}
# same code, with gCrown.y replace Live.y
{
  fit2 <- glm(Live.y~ldbh + as.factor(fsLevel) + Species.x,data=t12s,family='binomial')
  fit1 <- glm(Live.y~ldbh + as.factor(fsLevel) + Species.x+ as.factor(fsLevel):Species.x,data=t12s,family='binomial')
  BIC(fit1)
  BIC(fit2)
  summary(fit1)
  summary(fit2)
  
  nvals <- 11
  ldbh.vals <- seq(-0.5,2,length.out=nvals)
  fsLevels.vals <- c(0:3)
  nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.x=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  spname='QUEAGR'
  plotSpecies <- function(spname,nd.tmp=nd2) {
    tmp <- t12s[which(t12s$Species.x==spname),]
    #op=par(mfrow=c(1,2))
    
    # plot(tmp$FireSev,tmp$Live.y,main=spname,xlim=range(t12s$FireSev,na.rm=T))
    # i=1
    # for (i in 1:nvals) {
    #   ndt <- nd2[which(nd2$Species.x==spname & nd2$ldbh==ldbh.vals[i]),]
    #   lines(ndt$FireSev,ndt$pSurvAll)
    #   text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
    # }
    
    plot(tmp$ldbh,tmp$Live.y,xlim=range(t12s$ldbh,na.rm=T))
    i=2
    for (i in 1:length(unique(nd.tmp$fsLevel))) {
      ndt <- nd.tmp[which(nd2$Species.x==spname & nd.tmp$fsLevel==fsLevels.vals[i]),]
      lines(ndt$ldbh,ndt$pSurvAll)
      text(ldbh.vals[5],ndt$pSurvAll[which(ndt$ldbh==ldbh.vals[5])],fsLevels.vals[i])
    }
    #par(op)
  }
  plotSpecies('QUEGAR')
  
  # obtain predicted main effects of fire severity for each species at a common size
  nvals <- 1
  ldbh.vals <- log10(2)
  fsLevels.vals <- c(0:3)
  nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,length(fsLevels.vals)),fsLevels=rep(fsLevels.vals,each=nvals)))
  dim(nd)
  head(nd)
  
  nd2 <- data.frame(Species.x=rep(spA,each=nrow(nd)),ldbh=rep(nd$ldbh,length(spA)),fsLevel=rep(nd$fsLevels,length(spA)))
  
  nd2$pSurvAll <- predict(fit2,newdata=nd2,type='response')
  head(nd2)
  
  barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==0),],ylim=c(0,1))
  barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==1),],ylim=c(0,1))
  barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==2),],ylim=c(0,1))
  barplot(pSurvAll ~ Species.x,data=nd2[which(nd2$fsLevel==3),],ylim=c(0,1))
}
## END RECODE HERE

table(t12$Dead.y)
table(t12$Live.y)
table(t12$TB)
table(t12$gCrown.y)
table(t12$Live.y,t12$TB)
table(t12$Live.y,t12$gCrown.y)

## ANALYSIS FOR ONE SPECIES
spA

# pick one of these!
t12sp <- t12[which(t12$Species.x %in% c('QUEAGR')),]
t12sp <- t12[which(t12$Species.x %in% spA),]
t12sp <- t12[which(t12$Species.x %in% spA & t12$fsLevel>1),]

dim(t12sp)
t12sp <- t12sp[which(!is.na(t12sp$ldbh)),]
dim(t12sp)

# code to run a binomial model and plot response curve with data
# MORTALITY
nd <- with(t12sp,data.frame(ldbh=seq(min(t12sp$ldbh,na.rm=T),max(t12sp$ldbh,na.rm=T),length.out=101)))
nd$ldbh2 <- nd$ldbh^2

plot(t12sp$ldbh,t12sp$Dead.y)
fit <- glm(Dead.y~ldbh+ldbh2,data=t12sp,family='binomial')
nd$pMortality <- predict(fit,newdata=nd,type='response')
lines(nd$ldbh,nd$pMortality)

#TOPKILL WITH RESPROUT
plot(t12sp$ldbh,t12sp$TB)
fit <- glm(TB~ldbh+ldbh2 ,data=t12sp,family='binomial')
summary(fit)
nd$pTB <- predict(fit,newdata=nd,type='response')
lines(nd$ldbh,nd$pTB)

#GREEN CROWN
plot(t12sp$ldbh,t12sp$gCrown.y)
fit <- glm(gCrown.y~ldbh+ldbh2,data=t12sp,family='binomial')
summary(fit)
nd$pGCrown <- predict(fit,newdata=nd,type='response')
lines(nd$ldbh,nd$pGCrown)

#plot all three
plot(t12sp$ldbh,t12sp$gCrown.y)
lines(nd$ldbh,nd$pGCrown,col='green')
lines(nd$ldbh,nd$pTB,col='red')
lines(nd$ldbh,nd$pMortality)
nd$pTOT <- apply(nd[,c('pGCrown','pTB','pMortality')],1,sum)
summary(nd$pTOT)

# NOW FIT MULTINOMIAL
t12sp$PFstatus <- (-1)
t12sp$PFstatus[which(t12sp$Live.y==0)] <- 0
t12sp$PFstatus[which(t12sp$TB==1)] <- 1
t12sp$PFstatus[which(t12sp$gCrown.y==1)] <- 2
table(t12sp$PFstatus)

require(nnet)

# MULTINOMIAL - QUADRATIC CAN BE ADDED HERE '+ldbh2' - changes results some
fit1 <- multinom(as.factor(PFstatus) ~ ldbh +ldbh2 ,data=t12sp)
fit1
head(round(fitted(fit1),2))
dim(fitted(fit1))

plot(t12sp$ldbh,t12sp$Live.y)
points(t12sp$ldbh,fitted(fit1)[,1],col='black')
points(t12sp$ldbh,fitted(fit1)[,2],col='red')
points(t12sp$ldbh,fitted(fit1)[,3],col='green')
summary(apply(fitted(fit1)[,1:3],1,sum))
#######################

# END CLEAN CODE HERE!!!


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
#fit <- glmer(Live.y~ldbh + FireSev + Species.x + (1 |Plot.x),data=t12s,family='binomial')

######## TOPKILL ANALYSIS
# current scoring has dead trees as topkilled. Change so topkill is only for those that are alive
t12$TopkillLive.y <- 0
t12$TopkillLive.y[which(t12$Topkill.y==1 & t12$Live.y==1)] <- 1

table(t12$Dead.y,t12$TopkillLive.y,t12$gCrown.y)
t12[which(t12$gCrown.y==0 & t12$TopkillLive.y==0 & t12$Dead.y==0),]
# TWO INDIVIDUALS WITH PROBLEM DATA: 3415, 4437 (already identified those above - if they've been fixed and don't show up at this point, delete this line)

t12a <- t12[which(t12$FireSev>100),]

dim(t12a)
names(t12a)
head(t12a)
t12a$ldbh2 <- t12a$ldbh^2

op=par(mfrow=c(1,1))
plot(Dead.y~ldbh,data=t12a)
fitD <- glm(Dead.y~ldbh+ldbh2,data=t12a,family='binomial')
fitD
nd <- data.frame(ldbh=seq(min(t12a$ldbh,na.rm=T),max(t12a$ldbh,na.rm=T),length.out=101))
nd$ldbh2 <- nd$ldbh^2
head(nd)
nd$pDead <- predict(fitD,nd,type='response')
lines(pDead~ldbh,data=nd,lwd=2,col='black')

#plot(Topkill.y~ldbh,data=t12a)
fitT <- glm(TopkillLive.y~ldbh+ldbh2,data=t12a,family='binomial')
fitT
nd$pTopKill <- predict(fitT,nd,type='response')
lines(pTopKill~ldbh,data=nd,lwd=2,col='red')

#plot(gCrown.y~ldbh,data=t12a)
fitG2 <- glm(gCrown.y~ldbh+ldbh2,data=t12a,family='binomial')
fitG2
BIC(fitG2)
fitG1 <- glm(gCrown.y~ldbh,data=t12a,family='binomial')
fitG1
BIC(fitG1)

nd$pGreen <- predict(fitG1,nd,type='response')
lines(pGreen~ldbh,data=nd,lwd=2,col='green')

par(op)

#### NOW TRY MULTINOMIAL### SKIP THIS SECTION-###
require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(t12[,c('ldbh','PFstatus','Species.x','FireSev')],1,allNotNA)
table(rComp)
t12a <- t12[rComp,]

fit1 <- multinom(PFstatus ~ ldbh,data=t12a)
fit1
head(round(fitted(fit1),2))

plot(t12a$ldbh,t12a$Live.y)
points(t12a$ldbh,fitted(fit1)[,1],col='red')
points(t12a$ldbh,fitted(fit1)[,2],col='gray')
points(t12a$ldbh,fitted(fit1)[,3],col='green')


##########
head(t12)
dim(t12)
t12$dupStatusCheck <- apply(t12[,c('Dead.y','TopkillLive.y','gCrown.y')],1,sum,na.rm=T)
table(t12$dupStatusCheck)
## same 2 as above

###end skip section ###

#### PROVISIONALLY ASSIGN TO THREE CLASSES
t12$PFstatus <- (-1)
t12$PFstatus[which(t12$Live.y==0)] <- 0
t12$PFstatus[which(t12$TopkillLive.y==1)] <- 1
t12$PFstatus[which(t12$Topkill.y==0 & t12$gCrown.y==1)] <- 2
table(t12$PFstatus)

# now assign all remaing NAs to dead - TEMP STEP
t12$PFstatus[which(t12$PFstatus==(-1))] <- 0
table(t12$PFstatus)

## NOW TRY MULTINOMIAL
# setup discrete FireSev
summary(t12$FireSev)
t12$dFS <- cut(t12$FireSev,breaks = c(-200,35,130,298,1000))
table(t12$dFS)

require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(t12[,c('ldbh','PFstatus','Species.x','FireSev')],1,allNotNA)
table(rComp)
t12a <- t12[rComp,]

fit1 <- multinom(PFstatus ~ ldbh,data=t12a)
fit1
head(round(fitted(fit1),2))

plot(t12a$ldbh,t12a$Live.y)
points(t12a$ldbh,fitted(fit1)[,1],col='red')
points(t12a$ldbh,fitted(fit1)[,2],col='gray')
points(t12a$ldbh,fitted(fit1)[,3],col='green')

#Size and Fire Sev 
fit2 <- multinom(PFstatus ~ ldbh + dFS,data=t12a)
fit2
BIC(fit2)
head(round(fitted(fit2),2))

plot(t12a$ldbh,t12a$Live.y)
points(t12a$ldbh,fitted(fit2)[,1],col='red')
points(t12a$ldbh,fitted(fit2)[,2],col='grey')
points(t12a$ldbh,fitted(fit2)[,3],col='green')

# Add Species
fit3 <- multinom(PFstatus ~ ldbh + dFS + Species.x,data=t12a)
fit3
BIC(fit3)
head(round(fitted(fit3),2))

plot(t12a$ldbh,t12a$Live.y)
points(t12a$ldbh,fitted(fit3)[,1],col='red')
points(t12a$ldbh,fitted(fit3)[,2],col='grey')
points(t12a$ldbh,fitted(fit3)[,3],col='green')

# Add Species * size interaction
fit3x <- multinom(PFstatus ~ ldbh + dFS + Species.x + dFS:ldbh,data=t12a)
fit3x
BIC(fit3x)
head(round(fitted(fit3x),2))

plot(t12a$ldbh,t12a$Live.y)
points(t12a$ldbh,fitted(fit3x)[,1],col='red')
points(t12a$ldbh,fitted(fit3x)[,2],col='black')
points(t12a$ldbh,fitted(fit3x)[,3],col='green')

# COMPARE BIC
BIC(fit1)
BIC(fit2)
BIC(fit3)
BIC(fit3x)

# plot for selected values
rsel <- which(t12a$Species.x=='PSEMEN' & as.numeric(t12a$dFS)>0)
plot(t12a$ldbh[rsel],t12a$PFstatus[rsel])
points(t12a$ldbh[rsel],fitted(fit3x)[rsel,1],col='red',pch=19)
points(t12a$ldbh[rsel],fitted(fit3x)[rsel,2],col='black',pch=19)
points(t12a$ldbh[rsel],fitted(fit3x)[rsel,3],col='green',pch=19)

## logit model and plot
d=t12
xcn='ldbh'
ycn='Live.y'
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
logit2Plot(t12,'ldbh','Live.y')

