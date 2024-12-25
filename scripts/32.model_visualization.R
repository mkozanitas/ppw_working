## Visualizing model results
rm(list=ls())
local.dir <- '/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/Fire_2017/Demography paper 2024/model_fitting/'
mod <- 'brm'
spName <- 'QUEAGR'
fit.type <- 'Poly'
dep.var <- 'Multinom' #'H_live','H_gCxLive'

(fname <- paste(local.dir,paste(mod,spName,fit.type,dep.var,'rds',sep='.'),sep=''))
fit <- readRDS(fname)

attributes(fit)
head(fit$data)
fsLevels <- sort(unique(fit$data$fsCat2))
d10.17r <- range(fit$data$d10.17)

dtemp <- seq(d10.17r[1],d10.17r[2],length.out=100)
pd <- data.frame(d10.17=rep(dtemp,length(fsLevels)),fsCat=rep(fsLevels,each=length(dtemp)))
pd$fsCat2 <- factor(pd$fsCat)
pd$Plot <- sort(unique(fit$data$Plot))[1]
pd$TreeNum <- min(fit$data$TreeNum,na.rm=T)

pd$pval <- predict(fit,newdata = pd,type='response',allow_new_levels = T)
head(pd)
names(pd)

fsLevels

op=par(mfrow=c(1,3))
for (i in fsLevels) {
  rsel <- which(pd$fsCat2==i)
  plot(pd$d10.17[rsel],pd$pval[rsel,1],ylim=range(pd$pval))
  points(pd$d10.17[rsel],pd$pval[rsel,2],col='red')
  points(pd$d10.17[rsel],pd$pval[rsel,3],col='green')
}
par(op)

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
