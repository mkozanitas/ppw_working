## Visualizing model results
rm(list=ls())
library(rstantools)
source('scripts/31.functionsForAnalysis.R')
unlockBinding("last.warning", baseenv())

mod <- 'brm'
fit.type <- c('MN.Quad','MN.Splk3','MN.Splk6','MN.Splk20')
f <- 2
dep.var <- 'fate3.18' #'Live.18','gCxLv', 'fate3.18', 'gCrown.18'
iter <- 'i50000'
uhn <- 'Hect' # or 'Plot', if no FS: 'Hect.noFS' or 'Plot.noFS'

# UMBCAL QUEGAR QUEAGR PSEMEN HETARB QUEKEL ARBMEN QUEDOU AMOCAL
spList <- c('UMBCAL','QUEGAR','QUEAGR','HETARB','QUEKEL','ARBMEN','QUEDOU','AMOCAL','EHRO','WHTO','R.Shrub')

makeSpTabl <- function(dd) {
  fsCatLevels <- c('0.U','12.LM','3.H')
  SizeClassLevels <- c('SA','ST','MT','LT')
  spTab <- data.frame(fsCat=rep(fsCatLevels,each=4),SizeClass=rep(SizeClassLevels,3),N=NA,N.LT=NA,DN=NA,DR=NA,GC=NA)
  fc <- fsCatLevels[2]
  sc <- SizeClassLevels[1]
  for (fc in fsCatLevels)
    for (sc in SizeClassLevels)
    {
      st.row <- which(spTab$fsCat==fc & spTab$SizeClass==sc)
      if (sc=='LT') {
        rsel <- which(dd$fac.fsCat==fc & dd$Type.17.4 == 'LT')
        spTab$N.LT[st.row] <- length(rsel)
        rsel <- which(dd$fac.fsCat==fc & dd$Type.17.4 %in% c('LT','HLT'))
      } else rsel <- which(dd$fac.fsCat==fc & dd$Type.17.4 == sc)
      ddt <- dd[rsel,]
      spTab$N[st.row] <- length(rsel)
      if (length(rsel)!=0) {
        spTab$DN[st.row] <- length(dd$fate3.18[which(ddt$fate3.18=='DN')])/spTab$N[st.row]
        spTab$DR[st.row] <- length(dd$fate3.18[which(ddt$fate3.18=='DR')])/spTab$N[st.row]
        spTab$GC[st.row] <- length(dd$fate3.18[which(ddt$fate3.18=='GC')])/spTab$N[st.row]
        spTab[st.row,c('DN','DR','GC')] <- round(spTab[st.row,c('DN','DR','GC')],3)
      }
    }
  return(spTab)
}

sink('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/Fire_2017/Demography paper 2024/model_figs_2025/species-model-output.txt')
spN <- 1
for (spN in 1:length(spList))
{
  print(spN)
  spName <- spList[spN]
  dd <- readRDS(paste(local.dir,'/',paste(mod,spName,'dd.rds',sep='.'),sep=''))
  dim(dd)
  names(dd)
  dd$large.trees <- 'Sap+ST'
  dd$large.trees[which(dd$d10.17>=1.39076)] <- 'LT'
  
  #reset.warnings()
  (mfname <- paste(local.dir,'/',paste(mod,spName,uhn,fit.type[f],dep.var,iter,'rds',sep='.'),sep=''))
  (wfname <- paste(local.dir,'/',paste(mod,spName,uhn,fit.type[f],dep.var,iter,'WARNINGS.rds',sep='.'),sep=''))
  #mf <- multifit1
  mf <- readRDS(mfname)
  wf <- readRDS(wfname)
  
  rstantools::bayes_R2(mf)
  
  print(c(spName,nrow(dd)))
  print(mf)
  print(wf)
  
  # sample size and fate
  #print(table(dd$fate3.18,dd$fac.fsCat))
  spTab <- makeSpTabl(dd)
  write.csv(spTab,paste('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/Fire_2017/Demography paper 2024/model_figs_2025/',spName,'spTable.csv',sep=''))
  
  #print(table(dd$Survey,dd$large.trees))
  xlims <- range(dd$d10.17,na.rm=T)
  
  #conditional_effects(mf, categorical=TRUE)
  # use xlims set here for local range, default is study wide range
  visualizeMultifitBayes(mf,sp=fit.type[f],print.to.pdf=T,xlims=xlims) 
  #visualizeMultifitBayes(mf,sp=fit.type[f],print.to.pdf=T) 
}


# NOW FOR PSEMEN
reset.warnings()
spName <- 'PSEMEN'
dd <- readRDS(paste(local.dir,'/',paste(mod,spName,'dd.rds',sep='.'),sep=''))
dim(dd)
names(dd)
dd$large.trees <- 'Sap+ST'
dd$large.trees[which(dd$d10.17>=1.39076)] <- 'LT'
table(dd$Survey)

(mfname <- paste(local.dir,'/brm.PSEMEN.TRUE.10000.BERN.Splk3.Live18.rds',sep=''))
(wfname <- paste(local.dir,'/brm.PSEMEN.TRUE.10000.BERN.Splk3.Live18.WARNINGS.rds',sep=''))

spTab <- makeSpTabl(dd)
write.csv(spTab,paste('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/Fire_2017/Demography paper 2024/model_figs_2025/',spName,'spTable.csv',sep=''))

print(c(spName,nrow(dd)))

mf <- fit5brm.noplot
mf <- readRDS(mfname)
wf <- readRDS(wfname)
print(mf)
print(wf)

xlims <- range(dd$d10.17,na.rm=T)

visualizeBernfitBayes(mf,sp='SPLK3',print.to.pdf=T,xlims=xlims) 
sink()

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
