## Visualizing model results
rm(list=ls())
source('scripts/31.functionsForAnalysis.R')

mod <- 'brm'
fit.type <- c('MN.Quad','MN.Splk3','MN.Splk6','MN.Splk20')
f <- 2
dep.var <- 'fate3.18' #'Live.18','gCxLv', 'fate3.18', 'gCrown.18'
iter <- 'i50000'

spList <- c('QUEAGR','UMBCAL','HETARB','AMOCAL','QUEGAR','ARBMEN','EHRO','WHTO','R.Shrub')
spName <- spList[6]
dd <- readRDS(paste(local.dir,'/',paste(mod,spName,'dd.rds',sep='.'),sep=''))

reset.warnings()
(mfname <- paste(local.dir,'/',paste(mod,spName,fit.type[f],dep.var,iter,'rds',sep='.'),sep=''))
(wfname <- paste(local.dir,'/',paste(mod,spName,fit.type[f],dep.var,iter,'WARNINGS.rds',sep='.'),sep=''))
mf <- readRDS(mfname)
wf <- readRDS(wfname)
print(wf)

visualizeMultifitBayes(mf,sp=fit.type[f]) 

# NOW FOR PSEMEN
reset.warnings()
spName <- 'PSEMEN'
(mfname <- paste(local.dir,'brm.PSEMEN.BERN.Splk3.Live18.rds',sep=''))
(wfname <- paste(local.dir,'brm.PSEMEN.BERN.Splk3.Live18.WARNINGS.rds',sep=''))
mf <- readRDS(mfname)
wf <- readRDS(wfname)
print(wf)

visualizeBernfitBayes(mf,sp='SPLK3') 

## code below just to review warnings
i=6
for (i in 1:length(spList))
{
  reset.warnings()
  f <- 2
  spName <- spList[i]
  (wfname <- paste(local.dir,paste(mod,spName,fit.type[f],dep.var,'WARNINGS.rds',sep='.'),sep=''))
  wf <- readRDS(wfname)
  print(spName)
  print(wf)
}
(wfname <- paste(local.dir,'brm.PSEMEN.BERN.Splk3.Live18.WARNINGS.rds',sep=''))
wf <- readRDS(wfname)
print('PSEMEN')
print(wf)
### END HERE FOR NOW

fit.type <- 'Poly' #'MN.Poly', 'Poly'
dep.var <- 'Live.18' #'Live.18','gCxLv', 'fate3.18', 'gCrown.18'
(fname <- paste(local.dir,paste(mod,spName,fit.type,dep.var,'rds',sep='.'),sep=''))
fit.Live <- readRDS(fname)

dep.var <- 'gCrown.18' #'Live.18','gCxLv', 'fate3.18', 'gCrown.18'
(fname <- paste(local.dir,paste(mod,spName,fit.type,dep.var,'rds',sep='.'),sep=''))
fit.gCrown <- readRDS(fname)

dep.var <- 'gCxLv' #'Live.18','gCxLv', 'fate3.18', 'gCrown.18'
(fname <- paste(local.dir,paste(mod,spName,fit.type,dep.var,'rds',sep='.'),sep=''))
fit.gCxLv <- readRDS(fname)

fit.type <- 'MN.Poly'
dep.var <- 'fate3.18' #'Live.18','gCxLv', 'fate3.18', 'gCrown.18'
(fname <- paste(local.dir,paste(mod,spName,fit.type,dep.var,'rds',sep='.'),sep=''))
fit.fate3 <- readRDS(fname)



fsLevels <- sort(unique(fit.Live$data$fac.fsCat))
dtemp <- seq(mnd,mxd,length.out=100)
pd <- data.frame(d10.17=rep(dtemp,length(fsLevels)),fsCat=rep(fsLevels,each=length(dtemp)))
pd$fac.fsCat <- factor(pd$fsCat)
pd$Plot <- sort(unique(fit.Live$data$Plot))[1]
pd$TreeNum <- min(fit.Live$data$TreeNum,na.rm=T)

pd$p.Live <- predict(fit.Live,newdata = pd,type='response',allow_new_levels = T)[,1]
pd$p.gCrown <- predict(fit.gCrown,newdata = pd,type='response',allow_new_levels = T)[,1]
pd$p.gCxLv <- predict(fit.gCxLv,newdata = pd,type='response',allow_new_levels = T)[,1]
pd$p.gCxLv.x.p.Live <- pd$p.gCxLv * pd$p.Live

plot(pd$p.gCrown,pd$p.gCxLv.x.p.Live)
abline(0,1)

pd$diff <- pd$p.gCxLv.x.p.Live - pd$p.gCrown
head(pd)

which.min(pd$diff)
pd[which.min(pd$diff),]
pd[which.max(pd$diff),]

plot(pd$p.Live,pd$diff);abline(h=0,lty=2)
plot(pd$fac.fsCat,pd$diff);abline(h=0,lty=2)
plot(pd$d10.17,pd$diff);abline(h=0,lty=2)

plot(pd$d10.17,pd$p.gCrown,ylim=c(0,1))
points(pd$d10.17,pd$p.Live,col='green')
points(pd$d10.17,pd$p.gCxLv.x.p.Live,col='red')

### BELOW NOT WORKING NOW



## now try plotting hierarchical model
dv1 <- 'H_live'
dv2 <- 'H_gCxLive'
(fname1 <- paste(local.dir,paste(mod,spName,fit.type,dv1,'rds',sep='.'),sep=''))
fit1 <- readRDS(fname1)
(fname2 <- paste(local.dir,paste(mod,spName,fit.type,dv2,'rds',sep='.'),sep=''))
fit2 <- readRDS(fname2)

names(pd)
pdh <- pd[,1:5]
dim(pdh)
pdh$pLive <- predict(fit1,newdata = pdh,type='response',allow_new_levels = T)
pdh$pgCxLive <- predict(fit2,newdata = pdh,type='response',allow_new_levels = T)
pdh$pDN <- 1-pdh$pLive
pdh$pGC <- pdh$pLive * pdh$pgCxLive
pdh$pDR <- 1-pdh$pDN-pdh$pGC

op=par(mfrow=c(1,3))
for (i in fsLevels) {
  rsel <- which(pdh$fsCat2==i)
  plot(pdh$d10.17[rsel],pdh$pDN[rsel],ylim=c(0,1))
  points(pdh$d10.17[rsel],pdh$pgCxLive[rsel],col='blue')
  points(pdh$d10.17[rsel],pdh$pLive[rsel],col='lightblue')  
  points(pdh$d10.17[rsel],pdh$pDR[rsel],col='red')
  points(pdh$d10.17[rsel],pdh$pGC[rsel],col='green')
}
par(op)
