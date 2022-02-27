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
allIndv <- readRDS('data/allIndv.Rdata')
head(allIndv)

## get percent survival by species and type, from 13-18
use.species <- spNames$x[-c(17,19,21,23,25,34,35,36)]
fst <- data.frame(Species=rep(use.species,each=2),Type=rep(c('SA','TR'),length(use.species)),N13=NA,N18.dead=NA,N18.TKB=NA,N18.BC=NA,N18.C=NA)
head(fst)                  

# select which individuals to eliminate - this can be done separately for each analysis
# take any individuals where plot and species match for the first two years
nrow(allIndv)
plot.ok <- allIndv$Num[which(allIndv$P13==allIndv$P18)] 
spec.ok <- allIndv$Num[which(allIndv$S13==allIndv$S18)] 
ps.ok <- intersect(plot.ok,spec.ok)
length(ps.ok)

# not missed in 2018 survey
nm <- all.id[[1]]$Num[which(all.id[[1]]$NA18==0)]
head(nm)
length(nm)

t1 <- all.id[[1]][which(intersect(ps.ok,nm) %in% all.id[[1]]$Num),]
t1 <- t1[order(t1$Num),]
nrow(t1)
t2 <- all.id[[2]][which(intersect(ps.ok,nm) %in% all.id[[2]]$Num),]
nrow(t2)
t2 <- t2[order(t1$Num),]

# UNMMERGED - ANALYSES BASED ON SHARED NUMS
i=3
for (i in 1:nrow(fst))
{
  sp <- fst$Species[i]
  ty <- fst$Type[i]
  init <- t1$Num[which(t1$Species==sp & t1$Type==ty & t1$Live==1)]
  fst$N13[i] <- length(init)
  
  ## If running line by line, look at the initial plants in t2. E.g. for i=3, why so much missing data?
  t2[t2$Num %in% init,]
  
  fin <- intersect(init,t2$Num[which(t2$Live==1)])
  # fin is the number of original ty surviving, whether or not they transitioned from SA->TR (or the other way!)
  fst$N18.dead[i] <- length(init) - length(fin)
  
  fin <- intersect(init,t2$Num[which(t2$bSprout==1 & t2$Topkill==1)])
  fst$N18.TKB[i] <- length(fin) 
  
  fin <- intersect(init,t2$Num[which(t2$bSprout==1 & t2$gCrown==1)])
  fst$N18.BC[i] <- length(fin) 
  
  fin <- intersect(init,t2$Num[which(t2$bSprout==0 & t2$gCrown==1)])
  fst$N18.C[i] <- length(fin) 
}
head(fst)
tail(fst)
fst <- fst[-which(fst$N13==0),]
fst

## count missing or double counted indviduals between t1 and t2
fst$nDiff <- fst$N13-apply(fst[,4:7],1,sum,na.rm=T)
fst[,c(1,3,8)]

fst$S13.18[which(fst$N13!=0)] <- fst$N18[which(fst$N13!=0)]/fst$N13[which(fst$N13!=0)]
fst[order(fst$N13,decreasing = T),]

sum(fst$N13)
sum(fst$N18)
sum(fst$N18)/sum(fst$N13)

fst
#barplot(t(as.matrix(fst[,c('N13','N18')])))

#### EXPAND ANALYSIS to ALL POST-FIRE FATES
head(t1)
table(t1$Live,t1$bSprout)

# NOW MERGE FOR ANALYSES BY INDIV
dim(t1)
dim(t2)

t12 <- merge(t1,t2,'Num',all = T)
dim(t12)
head(t12)
names(t12)
table(t12$Num)
head(which(is.na(t12$Plot.x)))

## choose fire severity metric
fsmet <- 'Tubbs.MTBS.RDNBR.30'
summary(fs[,fsmet])

f2t <- match(t12$Plot.x,fs$Plot)
head(f2t)
tail(f2t)
t12$FireSev <- fs[f2t,fsmet]
dim(t12)
tail(t12)

# CALC DBH - easier to connect with our field data and knowledge
t12$dbh <- ba2d(t12$Basal.Area.x)
summary(t12$dbh)

## USE LOG DBH for analysis
hist(t12$dbh)
t12$ldbh <- log10(t12$dbh)
hist(t12$ldbh)

#### SURVIVAL ANALYSIS

plot(t12$ldbh,t12$Live.y)

# **** mixing SA and TR - need to work on diameter conversion to get this right ***
N <- length(which(!is.na(t12$ldbh) & !is.na(t12$Live.y)))
fit <- glm(Live.y~ldbh,data=t12,family='binomial')
summary(fit)

nd <- with(t12,data.frame(ldbh=seq(-0.5,2,length.out=1001)))
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)

plot(t12$ldbh,t12$Live.y)
lines(nd$ldbh,nd$pSurvAll,lwd=4)
abline(h=0.5,lty=2)

#what is critical basal area to achieve 50% survival?
ld50 <- nd$ldbh[which(nd$pSurvAll>=0.5)[1]]
abline(v=ld50,lty=2)
(spRes <- data.frame(species='All',N=N,d50=ld50,slp=fit$coefficients[2]))

## now run by species for species with lots of data
spN <- table(t12$Species.x)
spN[order(spN)]
(spA <- names(spN)[which(spN>=25)])
spA <- spA[-which(spA=='QUEBER')] # drop QUEBER - none died!

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
spRes$d50 <- round(as.numeric(spRes$d50),3)
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
plotSP(t12,'QUEAGR')


## Full model with species
t12s <- t12[which(t12$Species.x %in% spA),]
(N <- length(which(!is.na(t12s$ldbh) & !is.na(t12s$Live.y))))
head(t12s)
head(fs)

fit <- glm(Live.y~ldbh * Species.x,data=t12s,family='binomial')
summary(fit)

nd <- with(t12s,data.frame(ldbh=seq(-0.5,2,length.out=101)))
nd <- with(t12s,data.frame(ldbh=0.5,Species.x=unique(Species.x)))
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)
barplot(nd$pSurvAll~nd$Species.x)

### MODELS WITH FIRE SEVERITY

fit <- glm(Live.y~ldbh + FireSev,data=t12s,family='binomial')
summary(fit)

nvals <- 11
ldbh.vals <- seq(-0.5,2,length.out=nvals)
FireSev.vals <- seq(min(fs[,fsmet],na.rm=T),max(fs[,fsmet],na.rm=T),length.out=nvals)
nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,nvals),FireSev=rep(FireSev.vals,each=nvals)))
head(nd)

nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)

plot(t12s$FireSev,t12s$Live.y)
i=1
for (i in 1:nvals) {
  ndt <- nd[which(nd$ldbh==ldbh.vals[i]),]
  lines(ndt$FireSev,ndt$pSurvAll)
  text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
}

plot(t12s$ldbh,t12s$Live.y)
i=1
for (i in 1:nvals) {
  ndt <- nd[which(nd$FireSev==FireSev.vals[i]),]
  lines(ndt$ldbh,ndt$pSurvAll)
  text(ldbh.vals[5],ndt$pSurvAll[which(ndt$ldbh==ldbh.vals[5])],FireSev.vals[i])
}

# model with size, fire, species
fit2 <- glm(Live.y~ldbh + FireSev + Species.x,data=t12s,family='binomial')
fit1 <- glm(Live.y~ldbh + FireSev + Species.x+ ldbh:Species.x,data=t12s,family='binomial')
BIC(fit1)
BIC(fit2)
summary(fit1)

nvals <- 11
ldbh.vals <- seq(-0.5,2,length.out=nvals)
FireSev.vals <- seq(min(fs[,fsmet],na.rm=T),max(fs[,fsmet],na.rm=T),length.out=nvals)
nd <- with(t12,data.frame(ldbh=rep(ldbh.vals,nvals),FireSev=rep(FireSev.vals,each=nvals)))
dim(nd)
head(nd)

nd2 <- data.frame(Species.x=rep(spA,each=121),ldbh=rep(nd$ldbh,11),FireSev=rep(nd$FireSev,11))
dim(nd2)

nd2$pSurvAll <- predict(fit1,newdata=nd2,type='response')
head(nd2)

plotSpecies <- function(spname) {
  tmp <- t12s[which(t12s$Species.x==spname),]
  #op=par(mfrow=c(1,2))
  
  plot(tmp$FireSev,tmp$Live.y,main=spname,xlim=range(t12s$FireSev,na.rm=T))
  i=1
  for (i in 1:nvals) {
    ndt <- nd2[which(nd2$Species.x==spname & nd2$ldbh==ldbh.vals[i]),]
    lines(ndt$FireSev,ndt$pSurvAll)
    text(FireSev.vals[5],ndt$pSurvAll[which(ndt$FireSev==FireSev.vals[5])],ldbh.vals[i])
  }
  
  plot(tmp$ldbh,tmp$Live.y,xlim=range(t12s$ldbh,na.rm=T))
  i=1
  for (i in 1:nvals) {
    ndt <- nd2[which(nd2$Species.x==spname & nd2$FireSev==FireSev.vals[i]),]
    lines(ndt$ldbh,ndt$pSurvAll)
    text(ldbh.vals[5],ndt$pSurvAll[which(ndt$ldbh==ldbh.vals[5])],FireSev.vals[i])
  }
  #par(op)
}

spA
png('figures/logit1.png',width = 800, height = 1200)
op=par(mfrow=c(4,2))
plotSpecies('AMOCAL')
plotSpecies('ARBMEN')
plotSpecies('ARCMAN')
plotSpecies('FRACAL')
par(op)
dev.off()

png('figures/logit2.png',width = 800, height = 1200)
op=par(mfrow=c(4,2))
plotSpecies('HETARB')
plotSpecies('PSEMEN')
plotSpecies('UMBCAL')
par(op)
dev.off()

png('figures/logit3.png',width = 800, height = 1200)
op=par(mfrow=c(4,2))
plotSpecies('QUEAGR')
plotSpecies('QUEDOU')
plotSpecies('QUEGAR')
plotSpecies('QUEKEL')
par(op)
dev.off()

### full model with random plot
#fit <- glmer(Live.y~ldbh + FireSev + Species.x + (1 |Plot.x),data=t12s,family='binomial')

######## TOPKILL ANALYSIS
plot(t12$ldbh,t12$Live.y)

t12a <- t12[which(t12$Live.y==1
                   & t12$FireSev>100),]
dim(t12a)
names(t12a)
head(t12a)
t12a$ldbh2 <- t12a$ldbh^2

op=par(mfrow=c(1,2))
plot(Topkill.y~ldbh,data=t12a)
fit <- glm(Topkill.y~ldbh+ldbh2,data=t12a,family='binomial')
fit
nd <- data.frame(ldbh=seq(min(t12a$ldbh,na.rm=T),max(t12a$ldbh,na.rm=T),length.out=101))
nd$ldbh2 <- nd$ldbh^2
head(nd)
nd$yPred <- predict(fit,nd,type='response')
lines(yPred~ldbh,data=nd)

plot(gCrown.y~ldbh,data=t12a)



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
