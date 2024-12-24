# File with functions sourced by '30.Analysis.for.ms'

d2ba <- function(x) pi*(x/2)^2

loadFireSeverity <- function(fsmet='Tubbs.MTBS.RDNBR.30',verbose=F)
{
  tail(fs)
  
  names(fs)
  fs$Plot
  fs$fsvar <- fs[,fsmet]
  if (verbose) print(summary(fs$fsvar))
  if (verbose) print(hist(fs$fsvar))
  
  # RDNBR fire severity levels
  # Unburned
  # Low <170, but burned
  # Medium 170-700
  # High: >700
  # unburned: 1308, 1309, 1311, 1312, 1327, 1344, 1347
  
  # Create a discrete FireSev variable
  hist(fs$fsvar,breaks=c(-100,0,50,100,150,200,300,400,500,600,700,800,1000))
  fs.breaks <- c(-100,135,430,1000) # Based on Parks et al. 2014, but not using intermediate split at 304
  # !!this is for RDNBR.30 !!
  fs$fsCat <- cut(fs$fsvar,fs.breaks)
  if (verbose) print(table(fs$fsCat))
  fs$fsLevel <- as.numeric(fs$fsCat)
  if (verbose) print(table(fs$fsLevel))
  
  # manually code unburned as 0
  unburned.plots <- c('PPW1308','PPW1309','PPW1311','PPW1312','PPW1327','PPW1344','PPW1347')
  fs$fsLevel[which(fs$Plot %in% unburned.plots)] <- 0
  if (verbose) print(table(fs$fsLevel))
  fs$fsCat <- as.factor(fs$fsLevel)
  if (verbose) print(table(fs$fsCat))
  if (verbose) print(str(fs$fsCat))
  return(fs)
  ####
}

convertHybrids <- function()
{
  # hard code conversion of QUEBEGA to QUEBER
  tAll$Species.13[which(tAll$Species.13=='QUEBEGA')] <- 'QUEBER'
  tAll$Species.17[which(tAll$Species.17=='QUEBEGA')] <- 'QUEBER'
  tAll$Species.18[which(tAll$Species.18=='QUEBEGA')] <- 'QUEBER'
  tAll$Species.19[which(tAll$Species.19=='QUEBEGA')] <- 'QUEBER'
  tAll$Species.20[which(tAll$Species.20=='QUEBEGA')] <- 'QUEBER'
  tAll$Species[which(tAll$Species=='QUEBEGA')] <- 'QUEBER'
  
  tAll$Species.13[which(tAll$Species.13=='QUEDOGA')] <- 'QUEDOU'
  tAll$Species.17[which(tAll$Species.17=='QUEDOGA')] <- 'QUEDOU'
  tAll$Species.18[which(tAll$Species.18=='QUEDOGA')] <- 'QUEDOU'
  tAll$Species.19[which(tAll$Species.19=='QUEDOGA')] <- 'QUEDOU'
  tAll$Species.20[which(tAll$Species.20=='QUEDOGA')] <- 'QUEDOU'
  tAll$Species[which(tAll$Species=='QUEDOGA')] <- 'QUEDOU'
  
  tAll$Species.13[which(tAll$Species.13=='QUEdec')] <- 'QUEDOU'
  tAll$Species.17[which(tAll$Species.17=='QUEdec')] <- 'QUEDOU'
  tAll$Species.18[which(tAll$Species.18=='QUEdec')] <- 'QUEDOU'
  tAll$Species.19[which(tAll$Species.19=='QUEdec')] <- 'QUEDOU'
  tAll$Species.20[which(tAll$Species.20=='QUEdec')] <- 'QUEDOU'
  tAll$Species[which(tAll$Species=='QUEdec')] <- 'QUEDOU'
  
  tAll <- tAll[-which(tAll$Species=='UNK'),]
  
  tAll$Type.13[which(tAll$Type.13=='A2')] <- 'SA'
  tAll$Type.17[which(tAll$Type.17=='A2')] <- 'SA'
  #1112.1
  tAll$Type.18[which(tAll$Num==1112.1)] <- 'SA'
  return(tAll)
}

drawTernaryPlots <- function()
{
  s2sg <- match(tAll$Species,spAtt$OrigSpecies)
  isna <- which(is.na(s2sg))
  tAll$Species[isna]
  tAll$SGroup <- spAtt$SGroup[s2sg]
  tAll$BA.17 <- d2ba(tAll$d10.17)
  
  SGBA <- tapply(tAll$BA.17,list(tAll$Plot,tAll$SGroup),sum,na.rm=T)
  SGBA[is.na(SGBA)] <- 0
  sumBA <- apply(SGBA,1,sum)
  pSGBA <- SGBA/sumBA
  pSGBA3 <- pSGBA
  pSGBA3[,4] <- pSGBA3[,3] + pSGBA3[,4]
  (pSGBA3 <- round(100*pSGBA3[,-3],1))
  
  for (i in c(F,T)) {
    if (i) pdf('results/ternary.pdf',8.8)
    TernaryPlot(alab='CON',blab='EHRO',clab='WHTO')
    TernaryPoints(pSGBA3,pch=19,col=fsCols[fs$fsLevel+1],cex=2)
    #TernaryText(pSGBA3,substr(row.names(pSGBA3),6,7)) 
    if (i) dev.off()
  }
  
  
  # classify plots
  pt <- data.frame(plot=row.names(pSGBA),vt=NA)
  pt$vt[which(pSGBA3[,1]>=50)] <- 'CON'
  pt$vt[which(pSGBA3[,2]>=50)] <- 'MH'
  pt$vt[which(pSGBA3[,3]>=50)] <- 'WO'
  #pt$vt[which(is.na(pt$vt) & pSGBA3[,1]<20)] <- 'MH-WO'
  #pt$vt[which(is.na(pt$vt) & pSGBA3[,3]<20)] <- 'MH-CON'
  pt$vt[which(is.na(pt$vt))] <- 'Mix3'
  pt$fsLevel <- fs$fsLevel
  pt
  write.csv(pt,'data/vegtypes-fs.csv')
  return(pt)
}

calcFatesTableBySpecies <- function()
{
  fst12 <- data.frame(SpCode=rep(use.species,each=2),Type=rep(c('SA','TR'),length(use.species)),N17=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
  head(fst12)                  
  tail(fst12)
  
  i=1
  for (i in 1:nrow(fst12))
  {
    sp <- fst12$SpCode[i]
    ty <- fst12$Type[i]
    temp <- tAll[which(tAll$Species==sp & tAll$Type.17==ty),]
    nrow(temp)
    
    fst12$N17[i] <- sum(temp$Live.17,na.rm=T)
    
    ## The next three lines are all equivalent - just using third one
    #fst12$N18.DN[i] <- length(which(temp$fate.18=='DN'))
    #fst12$N18.DN[i] <- length(which(temp$DN.18=='1'))
    fst12$N18.DN[i] <- sum(temp$DN.18,na.rm = T)
    
    fst12$N18.DR[i] <- sum(temp$DR.18,na.rm = T)
    fst12$N18.LN[i] <- sum(temp$LN.18,na.rm = T)
    fst12$N18.LR[i] <- sum(temp$LR.18,na.rm = T)
    miss <- which(temp$Live.17==1 & is.na(temp$DN.18)==1)
    if (length(miss)>0) for (j in 1:length(miss)) print(temp[miss[j],c('Plot','Num')])
    fst12$nMissing <- fst12$N17-(fst12$N18.DN+fst12$N18.DR+fst12$N18.LN+fst12$N18.LR)
  }
  return(fst12)
}

reduce_fst12 <- function()
{
  # reduce fst12 to common species
  fst12c <- fst12[which(fst12$SpCode %in% common.species),-8]
  fst12c
  
  # add other trees and other shrubs to fst12c to make table for paper?
  (oth.shrubs <- spAtt$Species[which(spAtt$Shrub.Tree=='S' & spAtt$Common=='No')])
  (oth.trees <- spAtt$Species[which(spAtt$Shrub.Tree=='T' & spAtt$Common=='No')])
  
  newrow <- apply(fst12[which(fst12$SpCode %in% oth.shrubs & fst12$Type=='SA'),3:7],2,sum)
  fst12c <- rbind(fst12c,c('OTHSHR','SA',newrow))
  
  newrow <- apply(fst12[which(fst12$SpCode %in% oth.shrubs & fst12$Type=='TR'),3:7],2,sum)
  fst12c <- rbind(fst12c,c('OTHSHR','TR',newrow))
  
  newrow <- apply(fst12[which(fst12$SpCode %in% oth.trees & fst12$Type=='SA'),3:7],2,sum)
  fst12c <- rbind(fst12c,c('OTHTRS','SA',newrow))
  
  newrow <- apply(fst12[which(fst12$SpCode %in% oth.trees & fst12$Type=='TR'),3:7],2,sum)
  fst12c <- rbind(fst12c,c('OTHTRS','TR',newrow))
  row.names(fst12c) <- 1:nrow(fst12c)
  return(fst12c)
}

barplotFates <- function(d=tAll,fs='all')
{
  if ('all' %in% fs) fslevels <- 'fs.all'
  if ('drop-high' %in% fs)
  {
    d$fsCat[which(d$fsCat==3)] <- 2
    fslevels <- 'fs.nohi'
    
  }
  if ('low-medium' %in% fs)
  {
    d$fsCat[which(d$fsCat==2)] <- 1
    fslevels <- 'fs.dm'
  }
  
  d$fsCat2 <- factor(d$fsCat)
  d$fsCat <- d$fsCat2
  table(d$fsCat)
  
  
  ffs <- table(d$fate.18,d$fsCat)
  ffs
  ft <- apply(ffs,2,sum)
  ffsp <- ffs
  for (i in 1:ncol(ffsp)) ffsp[,i] <- ffs[,i]/ft[i]
  ffsp
  if (nrow(ffsp)==4) 
  {
    ffsp[3,] <- ffsp[3,]+ffsp[4,]
    ffsp <- ffsp[-4,]
  }
  ffsp
  barplot(ffsp,beside=F)
  
  return(ffsp)
}

prepareForBarPlots <- function()
{
  tAllm <- tAll
  # Optional - remove saplings added in 2018, where we might be introducing detection bias towards small survivors
  newSap <- which(is.na(tAllm$Year.17) & tAllm$Year.18==2018 & tAllm$Type.18=='SA')
  length(newSap)
  tAllm$Num[newSap]
  table(tAllm$Plot[newSap])
  table(tAllm$fate.18[newSap],tAllm$fsCat[newSap])
  
  # OPTION: comment this in or out to exercise option
  # tAllm <- tAllm[-newSap,]
  
  # OPTION: if needed trim data by stem size
  # tAllm <- tAllm[which(tAllm$ld10>=c(0)),]
  
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
  #FireLevels <- c('Mod+High'); FVals <- 2:3
  #FireLevels <- c('Low'); FVals <- 1
  #FireLevels <- c('None'); FVals <- 0
  #FireLevels <- c('None:Low'); FVals <- 0:1 #changes FireSev range to all c(1:3)
  FireLevels <- c('All'); FVals <- 0:3 #changes FireSev range to all c(1:3)
  
  #### Subset to selected fire values and species
  tAllm <- tAllm[which(tAllm$fsCat %in% FVals),]
  dim(tAllm)
  
  # remove rows with no size data
  #tAllm <- tAllm[which(!is.na(tAllm$ld10)),]
  #dim(tAllm)
  
  # Translate type and size data to categories
  tAllm$TSizeCat <- NA
  tAllm$TSizeCat[which(tAllm$Type.17=='SA' & tAllm$d10.17 <= max(tAll$d10.17,na.rm=T))] <- 'SA'
  tAllm$TSizeCat[which(tAllm$Type.17=='TR' & tAllm$DBH_cm.17 <= max(tAllm$DBH_cm.17,na.rm=T))] <- 'TR3'
  tAllm$TSizeCat[which(tAllm$Type.17=='TR' & tAllm$DBH_cm.17 < 20)] <- 'TR2'
  tAllm$TSizeCat[which(tAllm$Type.17=='TR' & tAllm$DBH_cm.17 < 10)] <- 'TR1'
  table(tAllm$TSizeCat,useNA='always')
  
  tAllm$SSizeCat <- NA
  tAllm$SSizeCat[which(tAllm$Type.17=='SA')] <- 'SA'
  tAllm$SSizeCat[which(tAllm$Type.17=='TR')] <- 'TR'
  table(tAllm$SSizeCat,useNA='always')
  
  return(tAllm)
}

barplotNonSprouters <- function(print.to.pdf=c(F,T))
{
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  tdat <- list()
  
  for (p2p in print.to.pdf) {
    if (p2p) pdf(paste('results/fates-non-sprouters.pdf',sep=''),6,6)
    
    op=par(mfrow=c(2,2))
    spsel <- 'PSEMEN';spname <- spsel
    temp <- tAllm[which(tAllm$Species %in% spsel),]
    
    temp <- temp[which(temp$TSizeCat %in% c('SA','TR1','TR2','TR3')),]
    temp$SizeCat <- temp$TSizeCat
    (tot <- table(temp$SizeCat,temp$fsCat))
    
    barplot(tot,beside = T,col=tree.cols,main=paste(spname,'N'))
    #barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
    barplot(1-table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(spname,'Survival'),ylim=c(0,1))
    tdat[[1]] <- temp
    
    
    spsel <- 'ARCMAN';spname <- spsel
    temp <- tAllm[which(tAllm$Species %in% spsel),]
    
    temp <- temp[which(temp$SSizeCat %in% c('SA','TR')),]
    temp$SizeCat <- temp$SSizeCat
    
    # include the very few high severity in medium
    temp$fsCat[which(temp$fsCat==3)] <- 2
    (tot <- table(temp$SizeCat,temp$fsCat))
    
    barplot(tot,beside = T,col=shrub.cols,main=paste(spname,'N'))
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=shrub.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
    tdat[[2]] <- temp
    
    
    par(op)
    if (p2p) dev.off()
  }
  return(tdat)
}

barplotOneNonSprouter <- function(ss,print.to.pdf=c(F,T),skip.op=T,plot.live=T)
{
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  tree.names <- spAtt$Species[which(spAtt$Shrub.Tree=='T')]
  shrub.names <- spAtt$Species[which(spAtt$Shrub.Tree=='S')]
  
  d <- tAllm[which(tAllm$Species %in% ss),]
  dim(d)
  
  if (ss %in% tree.names) {
    for (p2p in print.to.pdf) {
      if (p2p) pdf(paste('results/fates-',ss,'.pdf',sep=''),10,3)
      if (skip.op) op=par(mfrow=c(1,4))
      d <- d[which(d$TSizeCat %in% c('SA','TR1','TR2','TR3')),]
      d$SizeCat <- d$TSizeCat
      nrow(d)
      
      (tot <- table(d$SizeCat,d$fsCat))
      fates <- table(d$SizeCat,d$fsCat,d$fate3.18)
      barplot(tot,beside=T,col=tree.cols,main=paste(ss,'Sample Sizes'))
      if (!plot.live) barplot(fates[,,1]/tot,beside = T,col=tree.cols,main=paste(ss,'Mortality'),ylim=c(0,1))
      barplot(0/tot,beside = T,col=tree.cols,main=paste(ss,'Topkill-Resprout'),ylim=c(0,1))
      barplot(fates[,,2]/tot,beside = T,col=tree.cols,main=paste(ss,'Green Crown'),ylim=c(0,1))
      if (plot.live) barplot(fates[,,2]/tot,beside = T,col=tree.cols,main=paste(ss,'Live'),ylim=c(0,1))
      if (skip.op) par(op)
    }
  } else {
    for (p2p in print.to.pdf) {
      if (p2p) pdf(paste('results/fates-',ss,'.pdf',sep=''),10,3)
      if (skip.op) op=par(mfrow=c(1,4))
      d <- d[which(d$SSizeCat %in% c('SA','TR')),]
      d$SizeCat <- d$SSizeCat
      nrow(d)
      
      (tot <- table(d$SizeCat,d$fsCat))
      barplot(tot,beside=T,col=tree.cols,main=paste(ss,'Sample Sizes'))
      barplot(table(d$SizeCat,d$fsCat,d$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(ss,'Mortality'),ylim=c(0,1))
      barplot(table(d$SizeCat,d$fsCat,d$fate3.18)[,,2]/tot,beside = T,col=tree.cols,main=paste(ss,'Topkill-Resprout'),ylim=c(0,1))
      barplot(table(d$SizeCat,d$fsCat,d$fate3.18)[,,3]/tot,beside = T,col=tree.cols,main=paste(ss,'Green-Crown'),ylim=c(0,1))
      if (skip.op) par(op)
    }
  }
  return(d)
}

barplotSprouterSpecies <- function(ss,print.to.pdf=c(F,T),skip.op=F,plot.live=F,ss.name=NA)
{
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  tree.names <- spAtt$Species[which(spAtt$Shrub.Tree=='T')]
  shrub.names <- spAtt$Species[which(spAtt$Shrub.Tree=='S')]
  if (is.na(ss.name)) ss.name <- ss
  
  d <- tAllm[which(tAllm$Species %in% ss),]
  if (ss[1] %in% c('QUEGAR','AMOCAL')) d$fsCat[which(d$fsCat==3)] <- 2
  dim(d)
  
  if (ss[1] %in% tree.names) {
    for (p2p in print.to.pdf) {
      if (p2p) pdf(paste('results/fates-',ss,'.pdf',sep=''),10,3)
      if (skip.op) op=par(mfrow=c(1,4))
      d <- d[which(d$TSizeCat %in% c('SA','TR1','TR2','TR3')),]
      d$SizeCat <- d$TSizeCat
      nrow(d)
      
      (tot <- table(d$SizeCat,d$fsCat))
      fates <- table(d$SizeCat,d$fsCat,d$fate3.18)
      barplot(tot,beside=T,col=tree.cols,main=paste(ss.name,'Sample Sizes'))
      if (!plot.live) barplot(fates[,,1]/tot,beside = T,col=tree.cols,main=paste(ss.name,'Mortality'),ylim=c(0,1))
      barplot(fates[,,2]/tot,beside = T,col=tree.cols,main=paste(ss.name,'Topkill-Resprout'),ylim=c(0,1))
      barplot(fates[,,3]/tot,beside = T,col=tree.cols,main=paste(ss.name,'Green-Crown'),ylim=c(0,1))
      if (plot.live) barplot(1-(fates[,,1]/tot),beside = T,col=tree.cols,main=paste(ss.name,'Live'),ylim=c(0,1))
      if (skip.op) par(op)
    }
  } else {
    for (p2p in print.to.pdf) {
      if (p2p) pdf(paste('results/fates-',ss,'.pdf',sep=''),10,3)
      if (skip.op) op=par(mfrow=c(1,4))
      d <- d[which(d$SSizeCat %in% c('SA','TR')),]
      d$SizeCat <- d$SSizeCat
      nrow(d)
      
      (tot <- table(d$SizeCat,d$fsCat))
      fates <- table(d$SizeCat,d$fsCat,d$fate3.18)
      barplot(tot,beside=T,col=shrub.cols,main=paste(ss.name,'Sample Sizes'))
      if (!plot.live) barplot(fates[,,1]/tot,beside = T,col=shrub.cols,main=paste(ss.name,'Mortality'),ylim=c(0,1))
      barplot(fates[,,2]/tot,beside = T,col=shrub.cols,main=paste(ss.name,'Topkill-Resprout'),ylim=c(0,1))
      barplot(fates[,,3]/tot,beside = T,col=shrub.cols,main=paste(ss.name,'Green-Crown'),ylim=c(0,1))
      if (plot.live) barplot(1-(fates[,,1]/tot),beside = T,col=shrub.cols,main=paste(ss.name,'Live'),ylim=c(0,1))
      if (skip.op) par(op)
    }
  }
  return(d)
}

barplotSprouters <- function(print.to.pdf=c(F,T))
{
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  tdat <- list()
  for (p2p in print.to.pdf) {
    tsets <- list()
    tset.names <- c()
    tsets[[1]] <- c('QUEDOU','QUEGAR')
    tset.names[1] <- 'White Oaks'
    tsets[[2]] <- c('ARBMEN','QUEAGR','QUEKEL','UMBCAL')
    tset.names[2] <- 'EHRO'
    
    if (p2p) pdf(paste('results/fates-sprouters.pdf',sep=''),10,8)
    op=par(mfrow=c(3,4))
    
    i=1
    for (i in 1:length(tsets))
    {
      spsel <- tsets[[i]]
      spname <- tset.names[i]
      temp <- tAllm[which(tAllm$Species %in% spsel),]
      nrow(temp)
      
      temp <- temp[which(temp$TSizeCat %in% c('SA','TR1','TR2','TR3')),]
      temp$SizeCat <- temp$TSizeCat
      nrow(temp)
      
      # move high severity sample into medium
      if (i==1) temp$fsCat[which(temp$fsCat==3)] <- 2
      (tot <- table(temp$SizeCat,temp$fsCat))
      
      tree.cols <- c('grey90','grey60','grey30','black')
      barplot(tot,beside=T,col=tree.cols,main=paste(spname,'Sample Sizes'))
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=tree.cols,main=paste(spname,'Topkill-Resprout'),ylim=c(0,1))
      barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,col=tree.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
      tdat[[i]] <- temp
    }
    
    tsets[[1]] <- c('AMOCAL','FRACAL','HETARB')
    tset.names[1] <- 'Rshrubs'
    
    spsel <- tsets[[1]]
    spname <- tset.names[1]
    temp <- tAllm[which(tAllm$Species %in% spsel),]
    
    temp <- temp[which(temp$SSizeCat %in% c('SA','TR')),]
    temp$SizeCat <- temp$SSizeCat
    nrow(temp)    
    
    # only two at high severity, move to medium for analysis
    temp$fsCat[which(temp$fsCat==3)] <- 2
    
    (tot <- table(temp$SizeCat,temp$fsCat))
    
    shrub.cols <- c('grey90','black')
    barplot(tot,beside=T,col=shrub.cols,main=paste(spname,'Sample Sizes'))
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=shrub.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,2]/tot,beside = T,col=shrub.cols,main=paste(spname,'Topkill-Resprout'),ylim=c(0,1))
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,3]/tot,beside = T,col=shrub.cols,main=paste(spname,'Green-Crown'),ylim=c(0,1))
    tdat[[3]] <- temp
    
    par(op)
    if (p2p) dev.off()
  }
  return(tdat)
}

fitFatesMod <- function(d,spName=NA,fs='all',logt=T,live.only=F)
{
  # fs=low-medium - combine low and medium to one level
  # fs=drop high - 0,1,2 only, and combine any 3s into 2
  # logt - use log diameter
  
  # model fitting
  P50r <- data.frame(Species=spName,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  #head(d)
  table(d$Species)
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  table(d$fsCat)
  if ('all' %in% fs) fslevels <- 'fs.all'
  if ('drop-high' %in% fs)
  {
    d$fsCat[which(d$fsCat==3)] <- 2
    fslevels <- 'fs.nohi'
    
  }
  if ('low-medium' %in% fs)
  {
    d$fsCat[which(d$fsCat==2)] <- 1
    fslevels <- 'fs.dm'
  }
  
  d$fsCat2 <- factor(d$fsCat)
  d$fsCat <- d$fsCat2
  table(d$fsCat)
  
  if (live.only) d <- d[-which(d$fate.18=='DN'),]
  
  d$yvar <- d[,yvar]
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yvar),]
  nrow(d)
  length(unique(d$TreeNum))
  #head(d)
  
  fit5 <- glmmTMB(yvar~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yvar~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yvar~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yvar~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (FALSE) 
  {
    fit5s <- brm(yvar ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    saveRDS(fit5s,paste('results/brm-',yvar,'-',spName,'-',fslevels,'-',logt,'.RDS',sep=''))
  } else {
    fit5s <- readRDS(paste('results/brm-',yvar,'-',spName,'-',fslevels,'-',logt,'.RDS',sep=''))
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F,yvar,logt=logt,drawverts=F)
  p50
  P50r[1,1:6] <- c(spName,yvar,p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,yvar,logt=logt,drawverts=F)
  p50
  P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitFates2StepsMod.brm <- function(d,spName=NA,fs='all',logt=T,live.only=F)
{
  # fs=low-medium - combine low and medium to one level
  # fs=drop high - 0,1,2 only, and combine any 3s into 2
  # logt - use log diameter
  # local.dir <- local file directory for storing models - too large for github
  local.dir <- '/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/Fire_2017/Demography paper 2024/model_fitting'
  
  # model fitting
  table(d$Species)
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  table(d$fsCat)
  if ('all' %in% fs) fslevels <- 'fs.all'
  if ('drop-high' %in% fs)
  {
    d$fsCat[which(d$fsCat==3)] <- 2
    fslevels <- 'fs.nohi'
    
  }
  if ('low-medium' %in% fs)
  {
    d$fsCat[which(d$fsCat==2)] <- 1
    fslevels <- 'fs.dm'
  }
  
  d$fsCat2 <- factor(d$fsCat)
  d$fsCat <- d$fsCat2
  table(d$fsCat)
  
  # fit multinomial first
  table(d$fate.18)
  d$fate3 <- d$fate.18
  d$fate3[which(d$fate3 %in% c('LN','LR'))] <- 'GC'
  dd <- d[complete.cases(d$fsCat2,d$d10.17,d$fate3),]
  dim(dd)
  table(dd$fate3)
  table(dd$fate3,dd$fsCat)
  
  multifit1 <- brm(fate3 ~ s(d10.17, k=20, by=fsCat2) + fsCat2 + (1|Plot) + (1|TreeNum), data=dd,
                   family="categorical", chains = 2,
                   cores = 2, 
                   seed=726, 
                   #backend="cmdstanr",
                   refresh=100,
                   control=list(adapt_delta=0.95))
  summary(multifit1)
  
  # File is too large for github, will require local storage. 
  saveRDS(fit5brm,paste(local.dir,'/brm.',spSel,'.','Poly.Multinom.rds',sep=''))
  
  #Now fit hierarchical polynomial, Live first
  yvar <- 'Live.18'
  d$yvar <- d[,yvar]
  dd <- d[complete.cases(d$fsCat2,d$d10.17,d$yvar,d$Plot,d$TreeNum),]
  table(dd$Plot)
  dim(dd)
  
  fit5brm <- brm(yvar ~ d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=dd, family= 'bernoulli')
  print(summary(fit5brm))
  
  # File is too large for github, will require local storage. 
  saveRDS(fit5brm,paste(local.dir,'/brm.',spSel,'.','Poly.H_Live.rds',sep=''))
  
  # trouble shooting convergence failure - is it due to random effects
  if (FALSE) {
    fit5brm1 <- brm(yvar ~ d10.17 * fsCat2 + (1|Plot), data=dd, family= 'bernoulli')
  }
  
  # Now analyze live only
  dd <- dd[-which(dd$fate.18=='DN'),]
  table(dd$gCrown.18)
  table(dd$Plot)

  yvar <- 'gCrown.18'
  dd$yvar <- dd[,yvar]
  table(dd$yvar,useNA='always')
  dd <- dd[complete.cases(d$fsCat2,d$d10.17,d$yvar,d$Plot,d$TreeNum),]
  dim(dd)
  
  fit5brm <- brm(yvar ~ d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family= 'bernoulli')
  print(summary(fit5brm))

  # File is too large for github, will require local storage. 
  saveRDS(fit5brm,paste(local.dir,'/brm.',spSel,'.','Poly.H_gCxLive.rds',sep=''))
  
  # trouble shooting convergence failure - is it due to random effects
  if (FALSE) {
    fit5brm1 <- brm( yvar ~ d10.17 * fsCat2 + (1|Plot), data=d, family= 'bernoulli')
  }
}

fitFates2StepsMod <- function(d,spName=NA,fs='all',logt=T,live.only=F)
{
  # fs=low-medium - combine low and medium to one level
  # fs=drop high - 0,1,2 only, and combine any 3s into 2
  # logt - use log diameter
  
  # model fitting
  table(d$Species)
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  table(d$fsCat)
  if ('all' %in% fs) fslevels <- 'fs.all'
  if ('drop-high' %in% fs)
  {
    d$fsCat[which(d$fsCat==3)] <- 2
    fslevels <- 'fs.nohi'
    
  }
  if ('low-medium' %in% fs)
  {
    d$fsCat[which(d$fsCat==2)] <- 1
    fslevels <- 'fs.dm'
  }
  
  d$fsCat2 <- factor(d$fsCat)
  d$fsCat <- d$fsCat2
  table(d$fsCat)
  
  yvar <- 'Live.18'
  d$yvar <- d[,yvar]
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yvar,d$Plot,d$TreeNum),]
  table(d$Plot)
  dim(d)
  
  #fit5 <- glmmTMB(yvar~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  #fit5 <- brm(yvar~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  fit5brm <- brm( yvar ~ d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family= 'bernoulli')
  print(summary(fit5))
  AIC(fit5)
  
  pr <- predict_response(fit5,terms=c('d10.17','fsCat2'),margin = 'empirical')
  plot(pr,limits=c(0,1))

  fit5_pred <- ggpredict(fit5, c("d10.17", "fsCat2"))
  fit5_plot <- plot(fit5_pred)
  fit5_plot
  
  #mod2 <- lm(mpg ~ wt  qsec  factor(gear), data = mtcars)
  plot_comparisons(fit5, variables = "d10.17", condition = c("fsCat2"),vcov=T,re.form=NA)
  
  # don't understand x-axis
  p50Live <- plotPredictedValuesFit5(fit5,F,yvar,logt=logt,drawverts=F,ret.pv=T)
  
  # trouble shooting convergence failure - is it due to random effects
  if (FALSE) {
    fit5xi <- glmmTMB(yvar~d10.17 * fsCat2 + (1|Plot), data=d, family='binomial')
    print(summary(fit5xi))
    plotPredictedValuesFit5(fit5xi,F,yvar,logt=logt,drawverts=F,ret.pv=F)
    
    fit5xp <- glmmTMB(yvar~d10.17 * fsCat2 + (1|TreeNum), data=d, family='binomial')
    print(summary(fit5xp))
    plotPredictedValuesFit5(fit5xp,F,yvar,logt=logt,drawverts=F,ret.pv=F)
  }
  
  # now analyze live only
  dd <- d[-which(d$fate.18=='DN'),]
  table(dd$gCrown.18)
  table(dd$Plot)
  dim(dd)
  yvar <- 'gCrown.18'
  dd$yvar <- dd[,yvar]
  table(dd$yvar,useNA='always')
  
  fit5 <- glmmTMB(yvar~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=dd, family='binomial',control = glmmTMBControl(optCtrl = list(iter.max = 1000, rel.tol = 1e-6)))
  #fit5 <- glmmTMB(yvar~d10.17 + fsCat2 , data=dd, family='binomial',control = glmmTMBControl(optCtrl = list(iter.max = 1000, rel.tol = 1e-6)))
  if (FALSE) {
    require(car)
    vif(fit5)
    plot(d10.17~fsCat2,data=dd)
  }
  print(summary(fit5))
  AIC(fit5)
  pr <- predict_response(fit5,terms=c('d10.17','fsCat2'),margin = 'empirical')
  plot(pr,)
  p50GC <- plotPredictedValuesFit5(fit5,F,yvar,logt=logt,drawverts=F,ret.pv=T)
  
  # trouble shooting convergence failure - is it due to random effects
  if (FALSE) {
    fit5xi <- glmmTMB(yvar~d10.17 * fsCat2 + (1|Plot), data=dd, family='binomial',control = glmmTMBControl(optCtrl = list(iter.max = 1000, rel.tol = 1e-6)))
    print(summary(fit5xi))
    plotPredictedValuesFit5(fit5xi,F,yvar,logt=logt,drawverts=F,ret.pv=F)
    
    fit5xp <- glmmTMB(yvar~d10.17 * fsCat2 + (1|TreeNum), data=dd, family='binomial',,control = glmmTMBControl(optCtrl = list(iter.max = 1000, rel.tol = 1e-6)))
    print(summary(fit5xp))
    plotPredictedValuesFit5(fit5xp,F,yvar,logt=logt,drawverts=F,ret.pv=F)
  }
  
  p50L <- p50Live[[2]]
  p50G <- p50GC[[2]]
  all(p50L$d10.17==p50G$d10.17)
  p50J <- p50L[,c(1,2,6)]
  names(p50J)[3] <- 'pLive'  
  
  p50J$pGCxLive <- p50G$pval
  head(p50J)
  p50J$pGC <- p50J$pLive * p50J$pGCxLive
  p50J$pDR <- p50J$pLive * (1-p50J$pGCxLive)
  p50J$pDN <- 1 - p50J$pLive
  head(p50J)
  
  op=par(mfrow=c(1,3))
  for (i in c(0,1,3))
  {
    d <- p50J[which(p50J$fsCat==i),]
    plot(pDN~d10.17,data=d,ylim=c(0,1),main=paste(spName,'FS =',i))
    points(pDR~d10.17,data=d,col='brown')
    points(pGC~d10.17,data=d,col='green')
  }
  par(op)
  # fit5_for_lrt <- glmmTMB(yvar~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  # summary(fit5_for_lrt)
  # print('fire severity effect')
  # print(anova(fit5, fit5_for_lrt, "LRT"))
  # 
  # fit5_for_lrt2 <- glmmTMB(yvar~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  # summary(fit5_for_lrt2)
  # print('size effect')
  # print(anova(fit5, fit5_for_lrt2, "LRT"))
  # 
  # fit5_for_lrt3 <- glmmTMB(yvar~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  # summary(fit5_for_lrt3)
  # print('severity*size interaction effect')
  # print(anova(fit5, fit5_for_lrt3, "LRT"))
  # 
  # if (FALSE) 
  # {
  #   fit5s <- brm(yvar ~ 
  #                  s(d10.17, by=fsCat2, k=20)+
  #                  fsCat2+
  #                  (1|Plot)+
  #                  (1|TreeNum),
  #                data=d,
  #                family="Bernoulli",
  #                chains = 2, cores = 2, seed=237, 
  #                #backend="cmdstanr",
  #                control=list(adapt_delta=0.99))
  #   conditional_effects(fit5s)
  #   saveRDS(fit5s,paste('results/brm-',yvar,'-',spName,'-',fslevels,'-',logt,'.RDS',sep=''))
  # } else {
  #   fit5s <- readRDS(paste('results/brm-',yvar,'-',spName,'-',fslevels,'-',logt,'.RDS',sep=''))
  # }
  # 
  # # second argument is F for glmmTMB models and T for brm models
  # p50 <- plotPredictedValuesFit5(fit5,F,yvar,logt=logt,drawverts=F)
  # p50
  # P50r[1,1:6] <- c(spName,yvar,p50)
  # p50 <- plotPredictedValuesFit5(fit5s,T,yvar,logt=logt,drawverts=F)
  # p50
  # P50r[1,7:10]<- p50
  # 
  # return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitMultiNomMod <- function()
{
  # fs=low-medium - combine low and medium to one level
  # fs=drop-high - 0,1,2 only, and combine any 3s into 2
  # logt - use log diameter
  
  head(d)
  table(d$Species)
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  table(d$fsCat)
  if ('all' %in% fs) fslevels <- 'fs.all'
  if ('drop-high' %in% fs)
  {
    d$fsCat[which(d$fsCat==3)] <- 2
    fslevels <- 'fs.nohi'
    
  }
  if ('low-medium' %in% fs)
  {
    d$fsCat[which(d$fsCat==2)] <- 1
    fslevels <- 'fs.dm'
  }
  d$fsCat2 <- factor(d$fsCat)
  d$fsCat <- d$fsCat2
  table(d$fsCat)
  
  table(d$fate.18)
  d$fate3 <- d$fate.18
  d$fate3[which(d$fate3 %in% c('LN','LR'))] <- 'GC'
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$fate3),]
  dim(d)
  table(d$fate3)
  table(d$fate3,d$fsCat)
  
  multifit1 <- brm(fate3 ~ s(d10.17, k=20, by=fsCat2) + fsCat2 + (1|Plot) + (1|TreeNum), data=d,
                   family="categorical", chains = 2, cores = 2, 
                   seed=726, 
                   #backend="cmdstanr",
                   refresh=100,
                   control=list(adapt_delta=0.95))
  summary(multifit1)
  
  summary(d$d10.17)
  summary(d$TreeNum)
  dtemp <- seq(0,max(d$d10.17,na.rm=T),length.out=1000)
  (dfsCat <- sort(unique(d$fsCat)))
  pd <- data.frame(d10.17=rep(dtemp,length(dfsCat)),fsCat=rep(dfsCat,each=length(dtemp)))
  pd$fsCat2 <- factor(pd$fsCat)
  pd$Plot <- sort(unique(d$Plot))[1]
  pd$TreeNum <- min(d$TreeNum,na.rm=T)
  
  pd$pval <- predict(multifit1,newdata = pd,type='response',allow_new_levels = T)
  head(pd)
  names(pd)
  
  stacked <- F
  fslevel <- 1
  fsr <- which(pd$fsCat==fslevel)
  if (stacked) yval <- apply(pd$pval[fsr,1:3],1,sum) else yval <- pd$pval[fsr,1]
  plot(pd$d10.17[fsr],yval,main=paste(spName,'- FS',fslevel),col='black',pch=19,
       ,ylim=c(0,1),xlim=c(0,2.2),xaxt='n',xlab='Diameter at Breast Height',ylab='Proportion (predicted value)',cex.lab=1.2)
  dbh.ax <- c(0.3,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150)
  ax.marks <- log10(dbh.ax*1.176 + 1.07)
  axis(1,at=ax.marks,labels=c('.....Saplings.....',c('1','2',NA,'4',NA,'6',NA,NA,NA,'10','20',NA,'40',NA,'60',NA,'80',NA,'100','150')))
  if (!spName=='PSEMEN') 
    if (stacked) 
    {
      points(pd$d10.17[fsr],pd$pval[fsr,2]+pd$pval[fsr,3],col='brown',pch=19)
      points(pd$d10.17[fsr],pd$pval[fsr,3],col='green',pch=19)
    } else {
      points(pd$d10.17[fsr],pd$pval[fsr,2],col='brown',pch=19)
      points(pd$d10.17[fsr],pd$pval[fsr,3],col='green',pch=19)
    }
}

plotPredictedValuesFit5 <- function(mod=fit5,pred.brm=F,yvar,logt=F,drawverts=t,ret.pv=F)
{
  summary(d$d10.17)
  summary(d$TreeNum)
  p50 <- c()
  dtemp <- seq(0,max(d$d10.17,na.rm=T),length.out=1000)
  dfsCat <- sort(unique(d$fsCat))
  pd <- data.frame(d10.17=rep(dtemp,length(dfsCat)),fsCat=rep(dfsCat,each=length(dtemp)))
  pd$fsCat2 <- factor(pd$fsCat)
  pd$Plot <- sort(unique(d$Plot))[1]
  #pd$Plot <- sample(unique(d$Plot),nrow(pd),replace=T)
  #pd$TreeNum <- min(d$TreeNum,na.rm=T)
  pd$TreeNum <- sample(unique(d$TreeNum),nrow(pd),replace=T)
  if (pred.brm) pd$pval <- predict(mod,newdata = pd,type='response',allow_new_levels = T) else pd$pval <- predict(mod,newdata = pd,type='response',allow.new.levels = T,re.form=NA)
  head(pd)
  names(pd)
  
  fsCols <- c('blue','brown','orange','red')
  if (pred.brm) pd$yval <- pd$pval[,1] else pd$yval <- pd$pval
  plot(pd$yval~pd$d10.17,type='n',ylim=c(0,1),main=paste(spName,yvar))
  i=1
  for (i in 1:4) {
    pdt <- pd[which(as.numeric(pd$fsCat)==(i-1)),]
    points(pdt$yval~pdt$d10.17,col=fsCols[i],type='b')
    p50[i] <- pdt$d10.17[which(pdt$yval>=0.5)[1]]
  }
  if (pred.brm) fpre <- 'brm' else fpre <- 'ml'
  fname <- paste('results/pval-',fpre,'-',yvar,'-',spName,'-',fslevels,'-',logt,'.RDS',sep='')
  saveRDS(pd,fname)
  if (drawverts) if (logt) abline(v=log10(c(0,2.25,12.8,25)),lty=2) else abline(v=c(0,2.25,12.8,25),lty=2)
  if (ret.pv) return(list(p50,pd)) else return(p50)
}

plotPredictedValuesByFS <- function(fpre='brm',spSel='UMBCAL') #fpre='ml'
{
  yvar <- 'gCrown.18'
  fname <- paste('results/pval-',fpre,'-',yvar,'-',spSel,'-',fslevels,'-',logt,'.RDS',sep='')
  gcpv <- readRDS(fname)
  if (!spSel=='PSEMEN')
  {
    yvar <- 'DR.18'
    fname <- paste('results/pval-',fpre,'-',yvar,'-',spSel,'-',fslevels,'-',logt,'.RDS',sep='')
    drpv <- readRDS(fname)
    yvar <- 'Live.18'
    fname <- paste('results/pval-',fpre,'-',yvar,'-',spSel,'-',fslevels,'-',logt,'.RDS',sep='')
    lvpv <- readRDS(fname)
    head(gcpv)
    head(drpv)
    head(lvpv)
    # check alignment of sizes
    all(gcpv$d10.17==drpv$d10.17)
    all(gcpv$d10.17==lvpv$d10.17)
    
    apv <- data.frame(d10.17=gcpv$d10.17,fsCat=gcpv$fsCat,gcval=gcpv$yval,drval=drpv$yval,lvval=lvpv$yval,gdval=gcpv$yval+drpv$yval)
  } else {
    apv <- data.frame(d10.17=gcpv$d10.17,fsCat=gcpv$fsCat,gcval=gcpv$yval,drval=NA,lvval=gcpv$yval,gdval=NA)
  }
  head(apv)
  
  stacked <- FALSE
  fslevel <- 3
  fsr <- which(apv$fsCat==fslevel)
  plot(apv$d10.17[fsr], apv$gcval[fsr],main=paste(spSel,'- FS',fslevel),col='green',pch=19,
       ,ylim=c(0,1),xlim=c(0,2.2),xaxt='n',xlab='Diameter at Breast Height',ylab='Proportion (predicted value)',cex.lab=1.2)
  dbh.ax <- c(0.3,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150)
  ax.marks <- log10(dbh.ax*1.176 + 1.07)
  axis(1,at=ax.marks,labels=c('.....Saplings.....',c('1','2',NA,'4',NA,'6',NA,NA,NA,'10','20',NA,'40',NA,'60',NA,'80',NA,'100','150')))
  if (!spSel=='PSEMEN') 
  {
    if (stacked) 
    {
      points(apv$d10.17[fsr],apv$lvval[fsr],col='brown',pch=19)
    } else {
      points(apv$d10.17[fsr],apv$drval[fsr],col='brown',pch=19)
      points(apv$d10.17[fsr],1-apv$lvval[fsr],col='black',pch=19)
      #points(apv$d10.17[fsr],apv$lvval[fsr]-apv$gcval[fsr],col='brown',pch=19)
    }
  } else {
    if (stacked)
    {
      
    } else {
      points(apv$d10.17[fsr],1-apv$lvval[fsr],col='black',pch=19)
    }
  }
}

fitModelsDiscreteFates <- function()
{
  # white oaks
  {
    temp <- tdat[[1]]
    fit0 <- multinom(as.factor(fate3.18) ~ fsCat * SizeCat * Species,data=temp)
    AIC(fit0)
    fit1 <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + fsCat:Species + SizeCat:Species,data=temp)
    AIC(fit1)
    fit2a <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + fsCat:Species,data=temp)
    AIC(fit2a)
    fit2b <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + SizeCat:Species,data=temp)
    AIC(fit2b)
    fit2c <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:Species + SizeCat:Species,data=temp)
    AIC(fit2c)
    fit3 <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species ,data=temp)
    AIC(fit3)
    fit4a <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat  ,data=temp)
    AIC(fit4a)
    fit4b <- multinom(as.factor(fate3.18) ~ fsCat + Species  ,data=temp)
    AIC(fit4b)
    fit4c <- multinom(as.factor(fate3.18) ~ SizeCat + Species  ,data=temp)
    AIC(fit4c)
    fit5 <- multinom(as.factor(fate3.18) ~ SizeCat,data=temp)
    AIC(fit5)
  }
  
  # EHRO
  {
    temp <- tdat[[2]]
    multifit1 <- brm(factor(fate3.18) ~ d10.17*mo(as.integer(fsCat))*Species, data=temp,
                     family="categorical", chains = 2, cores = 2, seed=726, 
                     #backend="cmdstanr",
                     refresh=100)
    
    fit0 <- multinom(as.factor(fate3.18) ~ fsCat * SizeCat * Species,data=temp)
    AIC(fit0)
    fit1 <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + fsCat:Species + SizeCat:Species,data=temp)
    AIC(fit1) # drop three way interaction
    
    # try all combinations of 2 way interactions
    fit2a <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + fsCat:Species,data=temp)
    fit2b <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + SizeCat:Species,data=temp)
    fit2c <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:Species + SizeCat:Species,data=temp)
    fit2d <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + SizeCat:Species,data=temp)
    fit2e <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:Species,data=temp)
    fit2f <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat,data=temp)
    
    AIC(fit2a)
    AIC(fit2b)
    AIC(fit2c) # BEST
    AIC(fit2d)
    AIC(fit2e)
    AIC(fit2f)
  }
  # best model - I don't think main effects should be dropped when terms are present in interactions
  fit2c <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:Species + SizeCat:Species,data=temp)
  
  # resprouting shrubs
  {
    temp <- tdat[[3]]
    fit3 <- multinom(as.factor(fate3.18) ~ fsCat * SizeCat * Species,data=temp)
    AIC(fit3)
    fit2 <- multinom(as.factor(fate3.18) ~ fsCat + SizeCat + Species + fsCat:SizeCat + fsCat:Species + SizeCat:Species,data=temp)
    AIC(fit2) # keep three way interaction, so I'll keep everything else
  }
  
}

delayedMortality <- function()
{
  table(tAll$fate.18)
  table(tAll$fate.19)
  table(tAll$fate.18,tAll$fate.19,useNA = 'always')
  # Came back to life!! Should these be reclassified in 2018 based on 2019
  # DN-DR = 78
  # DN-LN = 46
  
  print('file ressurect.csv written to data, with 124 plants recorded as live in 2019, but dead in 2018')
  resurrect <- which(tAll$fate.18=='DN' & (tAll$fate.19=='DR' | tAll$fate.19=='LN'))
  write.csv(tAll[resurrect,c('Plot','Num')],'data/resurrect.csv')
  
  # delayed mortality
  # DR-DN = 98
  # DR-NA = 17
  # LN-DN = 57
  # LN-NA = 3
  # LR-DN = 3
  dmort <- c(which(tAll$fate.18=='DR' & tAll$fate.19=='DN'),
             which(tAll$fate.18=='DR' & is.na(tAll$fate.19)),
             which(tAll$fate.18=='LN' & tAll$fate.19=='DN'),
             which(tAll$fate.18=='LN' & is.na(tAll$fate.19)),
             which(tAll$fate.18=='LR' & tAll$fate.19=='DN'))
  print(length(dmort))
  # Total = 178
  # out of how many living in 2018?
  print('proportions of delayed mortality, overall and by group')
  print('overall')
  print(length(dmort)/sum(tAll$Live.18,na.rm=T))
  
  dmort <- tAll[dmort,]
  
  print('by fire severity')
  print(table(dmort$fsCat))
  print(table(dmort$fsCat)/table(tAll$fsCat[which(tAll$Live.18==1)]))
  print('by type')
  print(table(dmort$Type.17))
  print(table(dmort$Type.17)/table(tAll$Type.17[which(tAll$Live.18==1)]))
  
  dms <- table(dmort$Species)
  dma <- table(tAll$Species[which(tAll$Live.18==1)])
  xx <- dms[match(names(dma),names(dms))]
  names(xx) <- names(dma)
  xx[which(is.na(xx))] <- 0
  print(sort(xx))
  yy <- round(xx/dma,3)
  yy <- cbind(yy,dma)
  yy <- yy[which(row.names(yy) %in% spAtt$Species[which(spAtt$Common=='Yes')]),]
  print('by species, common species only')
  print(yy[order(yy[,1]),])
}

basalDiameterHistograms <- function()
{
  cs <- spAtt$Species[which(spAtt$Common=='Yes')]
  tAllm$ld10.17 <- log10(tAllm$d10.17)
  
  dbks <- c(-1,-0.25,seq(0,1.5,0.125),2.5)
  hist(tAllm$ld10.17,breaks=dbks)
  
  ssi <- data.frame(Species=rep(cs,each=4),fsCat=rep(0:3,length(cs)))
  
  op=par(mfcol=c(4,1),mar=c(3,5,2,1))
  for (s in 1:length(cs))
  {
    sp <- cs[s]
    tmp <- tAllm[which(tAllm$Species==sp),]
    for (i in 0:3) {
      N <- length(which(tmp$fsCat==i))
      print(c(sp,i,N))
      if (N>5) {
        hs <- hist(tmp$ld10.17[tmp$fsCat==i],breaks=dbks,xlab='n',main=paste(sp,'fire level',i,'N:',N))
      } else plot(1:2,1:2,type='n',xlab='',ylab='',xaxt='n',yaxt='n',main=paste(sp,'fire level',i,'N:',N))
    }
  }
  par(mfcol=c(1,1),mar=c(5,5,3,3))
  
}

basalDiameterLogisticRegressions <- function()
{
  
  # Based on visual examination, chose which combinations work
  # AMOCAL, 0,1,2 (move 2 plants with fsCat=3 in to fsCat=2)
  # ARCMAN - too few small to parameterize, barplot tells the story
  # FRACAL 1 and 2
  # HETARB 0,1,2
  # ARBMEN all
  # PSEMEN all
  # QUEAGR all
  # QUEDOU too few small to parameterize, barplot tells the story
  # QUEGAR too few small to parameterize, barplot tells the story
  # QUEKEL no for 0, combine 1 and 2
  # UMBCAL all
  
  # run binomials for Live.18 and gCrown.18, and show predicted values for each fire sev and only the range of sizes observed, at each fire level!
  
  # remove QUEDOU anbd QUEGAR
  cs <- spAtt$Species[which(spAtt$Common=='Yes' & spAtt$Shrub.Tree=='T')]
  (cs <- cs[-which(cs %in% c('QUEDOU','QUEGAR'))])
  
  # critical diameters - basal = 1 cm, dbh = 1, 10, 20 - converted to basal
  # D10 = DBH * 1.176 + 1.070
  cd <- c(1,c(1,10,20,50)*1.176+1.07)
  fsLevels <- 0:3
  
  xx <- expand.grid(cd,fsLevels,cs)
  
  spRes <- data.frame(Species=xx[,3],fsCat=as.factor(xx[,2]),d10.17=xx[,1])
  spRes$cdr <- match(spRes$d10.17,cd)
  spRes$ld10.17 <- log10(spRes$d10.17)
  spRes$Live.18 <- NA
  spRes$gCrown.18 <- NA
  
  par(mfrow=c(2,1),mar=c(5,5,3,1))
  i=1
  for (i in 1:length(cs))
  {
    (sp <- cs[i])
    tmp <- tAllm[which(tAllm$Species==sp),]
    if (sp=='AMOCAL') {
      tmp$fsCat[which(tmp$fsCat==3)] <- 2
      fsLevels <- 0:2
    } else if (sp=='FRACAL') {
      tmp <- tmp[-which(tmp$fsCat==0),]
      fsLevels <- 1:2
    } else if (sp=='HETARB') {
      fsLevels=0:2
    } else if (sp=='QUEGAR') {
      tmp$fsCat[which(tmp$fsCat==3)] <- 2
      fsLevels=0:2
    } else if (sp=='QUEKEL') {
      tmp <- tmp[-which(tmp$fsCat==0),]
      tmp$fsCat[which(tmp$fsCat==2)] <- 1
      fsLevels=c(1,3)
      spRes <- spRes[-which(spRes$Species==sp & spRes$fsCat %in% c(0,2)),]
    } else fsLevels <- 0:3
    table(tmp$fsCat)
    
    y=2
    for (y in 1:2)
    {
      if (y==1) ynm <- 'Live.18' else ynm <- 'gCrown.18'
      tmp$yval <- tmp[,ynm]
      tmp$ld10.17 <- log10(tmp$d10.17)
      fit <- glm(yval~ld10.17+fsCat,data=tmp,family='binomial')  
      fit
      pdat <- data.frame(ld10.17=NA,fsCat=rep(fsLevels,each=100))
      f=0
      for (f in fsLevels) {
        sr <- range(tmp$ld10.17[which(tmp$fsCat==f)],na.rm=T)
        sims <- seq(sr[1],sr[2],length.out=100)
        pdat$ld10.17[which(pdat$fsCat==f)] <- sims
      }
      pdat$fsCat <- as.factor(pdat$fsCat)
      pdat$yval <- predict(fit,pdat,type='response')
      spRes[which(spRes==sp),ynm] <- round(predict(fit,spRes[which(spRes==sp),],type='response'),3)
      plot(yval~ld10.17,data=pdat,type='n',ylim=c(0,1),main=paste(sp,ynm),xlab='Basal Diameter (cm, log10)',ylab='')
      for (f in fsLevels) lines(yval~ld10.17,data=pdat[which(pdat$fsCat==f),],col=fsCols[f+1],lwd=2)
      abline(v=log10(cd)[c(1,2,3,5)],lty=2)
      if (sp=='QUEAGR') {
        #abline(h=0.5,lty=2)
        pdatqa <- pdat
      }
    }
  }
  par(mfcol=c(1,1),mar=c(5,5,3,3))
  
  # this shows critical thresholds for 50% crown survival in QUEAGR, but not using this in the paper
  if (FALSE) {
    head(pdatqa)
    scrit <- c()
    st <- which(pdatqa$yval>=0.5 & pdatqa$fsCat==1)
    scrit[1] <- pdatqa$ld10.17[st[1]]
    st <- which(pdatqa$yval>=0.5 & pdatqa$fsCat==2)
    scrit[2] <- pdatqa$ld10.17[st[1]]
    st <- which(pdatqa$yval>=0.5 & pdatqa$fsCat==3)
    scrit[3] <- pdatqa$ld10.17[st[1]]
    scrit
    10^scrit
    # D10 = DBH * 1.176 + 1.070
    # (D10-1.07)/1.176) - DBH
    (10^scrit-1.07)/1.176
  }
  
  # use this instead
  spRes[which(spRes$Species=='QUEAGR' & spRes$cdr==3),]
  spRes[which(spRes$Species=='PSEMEN' & spRes$cdr==3),]
  
  sprn <- data.frame(Species=cs,SpOrd <- c('1Am','0Pm','3Qa','4Qk','2Uc'))
  spRes$SpOrd <- sprn$SpOrd[match(spRes$Species,sprn$Species)]
  
  par(mfrow=c(4,2),mar=c(3,5,3,1))
  for (i in c(1,2,3,5)) {
    barplot(Live.18~fsCat+SpOrd,data=spRes[which(spRes$cdr==i),],beside=T,xlab='',ylab='Live',col=fsCols)
    barplot(gCrown.18~fsCat+SpOrd,data=spRes[which(spRes$cdr==i),],beside=T,xlab='',ylab='Green Crown',col=fsCols)
  }
  par(op.reset)
  
  par(mfrow=c(2,2),mar=c(3,5,3,1))
  spRes2 <- spRes[which(spRes$Species %in% c('PSEMEN','QUEAGR','QUEKEL')),]
  FS4 <- c('PSME','QUAG','QUKE')
  spRes2$FS4 <- FS4[match(spRes2$SpOrd,c('0Pm','3Qa','4Qk'))]
  for (i in c(1,2,3,5)) {
    barplot(Live.18~fsCat+FS4,data=spRes2[which(spRes2$cdr==i),],beside=T,xlab='',ylab='Live',col=fsCols,cex.lab=1.5,ylim=c(0,1))
  }
  par(op.reset)
}

basalResprouts <- function()
{
  names(tAllm)
  cs <- spAtt$Species[which(spAtt$Shrub.Tree %in% c('S','T') & spAtt$Common=='Yes' & spAtt$Resprout=='Y')]
  tAllc <- tAllm[which(tAllm$Species %in% cs & (tAllm$DR.18==1 | tAllm$LR.18==1)),]
  dim(tAllc)
  tAllc$Basal.Resprout.Height_cm.18 <- as.numeric(tAllc$Basal.Resprout.Height_cm.18)
  tAllc$Basal.Resprout.Height_cm.19 <- as.numeric(tAllc$Basal.Resprout.Height_cm.19)
  
  # fix two outlier - remove this after data file is fixed
  tAllc$Basal.Resprout.Height_cm.18[which.max(tAllc$Basal.Resprout.Height_cm.18)] <- NA
  tAllc$Basal.Resprout.Count.18[which.max(tAllc$Basal.Resprout.Count.18)] <- NA
  
  # change zeros to NAs to avoid including in analysis
  tAllc$Basal.Resprout.Height_cm.18[which(tAllc$Basal.Resprout.Height_cm.18==0)] <- NA
  tAllc$Basal.Resprout.Count.18[which(tAllc$Basal.Resprout.Count.18==0)] <- NA
  tAllc$Basal.Resprout.Height_cm.19[which(tAllc$Basal.Resprout.Height_cm.19==0)] <- NA
  tAllc$Basal.Resprout.Count.19[which(tAllc$Basal.Resprout.Count.19==0)] <- NA
  
  tAllc$lBasal.Resprout.Count.18 <- log10(tAllc$Basal.Resprout.Count.18+1)
  hist(tAllc$lBasal.Resprout.Count.18)
  hist(log10(tAllc$Basal.Resprout.Height_cm.18))
  plot(tAllc$lBasal.Resprout.Count.18,log10(tAllc$Basal.Resprout.Count.19+1))
  abline(0,1)
  
  fit <- lm(Basal.Resprout.Height_cm.18~Species*fsCat,data=tAllc)
  anova(fit)
  hist(fit$residuals)
  par(mfrow=c(2,1),mar=c(5,5,3,1))
  boxplot(Basal.Resprout.Height_cm.18~Species,data=tAllc,ylab='Resprout height (cm), 2018')
  boxplot(Basal.Resprout.Height_cm.18~fsCat,data=tAllc,xlab='Fire Severity',ylab='Resprout height (cm), 2018')
  
  length(which(tAllc$Basal.Resprout.Height_cm.18>0 & tAllc$fsCat==0))
  
  fit <- lm(lBasal.Resprout.Count.18~Species*fsCat,data=tAllc)
  anova(fit)
  hist(fit$residuals)
  par(mfrow=c(2,1),mar=c(5,5,3,1))
  boxplot(lBasal.Resprout.Count.18~Species,data=tAllc,ylab='Resprout count (+1,log10), 2018')
  boxplot(lBasal.Resprout.Count.18~fsCat,data=tAllc,xlab='Fire Severity',ylab='Resprout count (+1,log10), 2018')
  range(tAllc$Basal.Resprout.Count.18,na.rm=T)
  
  tAllc$RespGrowth <- tAllc$Basal.Resprout.Height_cm.19-tAllc$Basal.Resprout.Height_cm.18
  fit <- lm(RespGrowth~Species*fsCat,data=tAllc)
  anova(fit)
  hist(fit$residuals)
  boxplot(RespGrowth~Species,data=tAllc,ylab='Resprout growth, 2018-19 (cm)')
  abline(h=0,lty=2)
  boxplot(RespGrowth~fsCat,data=tAllc,xlab='Fire Severity',ylab='Resprout growth, 2018-19 (cm)')
  abline(h=0,lty=2)
  
  par(op.reset)
  plot(Basal.Resprout.Height_cm.19~Basal.Resprout.Height_cm.18,data=tAllc,xlab='Basal Resprought Height (cm), 2018',ylab='Basal Resprought Height (cm), 2019',asp=1,xlim=c(0,500),ylim=c(0,500),col=fsCols[as.numeric(tAllc$fsCat)],pch=19,cex=1.5)
  fit <- lm(Basal.Resprout.Height_cm.19~Basal.Resprout.Height_cm.18,data=tAllc)
  abline(fit)
  anova(fit)
  abline(0,1,lty=2)
  
  relGrowth <- tAllc$Basal.Resprout.Height_cm.19/tAllc$Basal.Resprout.Height_cm.18
  summary(relGrowth)
  
  # PSEMEN with basal resprouts - doesn't affect analyses above
  xx <- which(tAllm$Basal.Resprout.Height_cm.18>=0 & tAllm$Species=='PSEMEN')
  tAllm[xx,]
  tAllm$Basal.Resprout.Count.18[xx] <- NA
  tAllm$Basal.Resprout.Height_cm.18[xx] <- NA
  xx <- which(tAllm$Basal.Resprout.Height_cm.19>=0 & tAllm$Species=='PSEMEN')
  tAllm[xx,]
  tAllm$Basal.Resprout.Count.19[xx] <- NA
  tAllm$Basal.Resprout.Height_cm.19[xx] <- NA
}

new2019recruits <- function()
{
  r9 <- which(is.na(tAll$fate.18) & tAll$Type.19 %in% c('SA','TR') & tAll$Species!='BACPIL')
  print('number of new recruits, by type, species, fire severity')
  print(table(tAll$Type.19[r9]))
  print(table(tAll$Species[r9],tAll$Type.19[r9]))
  newr <- table(tAll$fsCat[r9],tAll$Type.19[r9])
  print(newr)
  s18 <- table(tAllm$fsCat[which(tAllm$Live.18==1)],tAllm$Type.19[which(tAllm$Live.18==1)])
  #newr/s18
  
  print('saplings by species and plot for unburned;low+medium;high severity')
  r9ubsa <- intersect(r9,which(tAll$fsCat %in% c(0) & tAll$Type.19=='SA'))
  print(table(tAll$Species[r9ubsa],tAll$Plot[r9ubsa]))
  
  r9lmssa <- intersect(r9,which(tAll$fsCat %in% c(1,2) & tAll$Type.19=='SA'))
  print(table(tAll$Species[r9lmssa],tAll$Plot[r9lmssa]))
  
  r9hssa <- intersect(r9,which(tAll$fsCat==3 & tAll$Type.19=='SA'))
  print(table(tAll$Species[r9hssa],tAll$Plot[r9hssa]))
  
  # What was the species composition of plots where Ceanothus germinated
  (cp <- unique(tAll$Plot[intersect(r9,which(tAll$Species=='CEACUN'))]))
  table(tAll$Species[which(tAll$Plot %in% cp[1] & tAll$Live.17==1)])
  table(tAll$Species[which(tAll$Plot %in% cp[1] & tAll$Live.18==1)])
  table(tAll$Species[which(tAll$Plot %in% cp[2] & tAll$Live.19==1)])
  table(tAll$Species[which(tAll$Plot %in% cp[2] & tAll$Live.17==1)])
  table(tAll$Species[which(tAll$Plot %in% cp[2] & tAll$Live.18==1)])
  table(tAll$Species[which(tAll$Plot %in% cp[2] & tAll$Live.19==1)])
  
  (cp <- unique(tAll$Plot[intersect(r9,which(tAll$Species=='AMOCAL'))]))
  unique(tAll$fsCat[tAll$Plot %in% c('PPW1348','PPW1350')])
  
  # tree size individuals recorded as new - some of these are too big to be new, and some have no dbh which seems strage - need to assess?
  r9t <- which(is.na(tAll$fate.18) & tAll$Type.19=='TR')
  length(r9t)
  r9tdbh <- tAll$DBH_cm.19[r9t]
}


# ATTEMPT AT GENERIC MODELS BUT IN THE END FIT THEM FOR EACH SPECIES
fitSpeciesModelsContSize.gCrown <-  function()
{
  fit5 <- glmmTMB(gCrown.18~d10.17*factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5 <- glmmTMB(gCrown.18~s(d10.17, by=fsCat2, k=20)*factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(gCrown.18~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print(anova(fit5, fit5_for_lrt, "LRT"))
  AIC(fit5_for_lrt)
  
  fit5_for_lrt2 <- glmmTMB(gCrown.18~factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  AIC(fit5_for_lrt2)
  
  fit5_for_lrt3 <- glmmTMB(gCrown.18~d10.17 + factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  AIC(fit5_for_lrt3)
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  summary(d$d10.17)
  summary(d$TreeNum)
  dtemp <- seq(min(d$d10.17,na.rm=T),max(d$d10.17,na.rm=T),length.out=100)
  dfsCat <- 0:3
  pd <- data.frame(d10.17=rep(dtemp,4),fsCat=rep(dfsCat,each=length(dtemp)))
  pd$Plot <- 'PPW1301'
  pd$TreeNum <- 1032
  pd$pval <- predict(fit5,newdata = pd,type='response',allow.new.levels = T)
  head(pd)
  
  fsCols <- c('blue','brown','orange','red')
  plot(pval~d10.17,data=pd,type='n',ylim=c(0,1),main='Green Crown')
  for (i in 1:4) {
    pdt <- pd[which(as.numeric(pd$fsCat)==(i-1)),]
    points(pval~d10.17,data=pdt,col=fsCols[i],type='b')
  }
  
}

fitSpeciesModelsContSize.Live <-  function()
{
  d$d10.17sq <- d$d10.17^2
  fit5 <- glmmTMB(Live.18~d10.17+factor(fsCat)+d10.17:factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5q <- glmmTMB(Live.18~d10.17sq+d10.17+factor(fsCat)+d10.17sq:factor(fsCat)+d10.17:factor(fsCat)+ (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5q))
  print(anova(fit5, fit5q, "LRT"))
  
  
  fit5_for_lrt <- glmmTMB(Live.18~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print(anova(fit5, fit5_for_lrt, "LRT"))
  AIC(fit5_for_lrt)
  
  fit5_for_lrt2 <- glmmTMB(Live.18~factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  AIC(fit5_for_lrt2)
  
  fit5_for_lrt3 <- glmmTMB(Live.18~d10.17 + factor(fsCat) + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  AIC(fit5_for_lrt3)
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  summary(d$d10.17)
  summary(d$TreeNum)
  dtemp <- seq(min(d$d10.17,na.rm=T),max(d$d10.17,na.rm=T),length.out=100)
  dfsCat <- 0:3
  pd <- data.frame(d10.17=rep(dtemp,4),fsCat=rep(dfsCat,each=length(dtemp)))
  pd$d10.17sq <- pd$d10.17^2
  pd$Plot <- 'PPW1301'
  pd$TreeNum <- 1032
  pd$pval <- predict(fit5,newdata = pd,type='response',allow.new.levels = T)
  head(pd)
  
  fsCols <- c('blue','brown','orange','red')
  plot(pval~d10.17,data=pd,type='n',ylim=c(0,1),main='Live')
  for (i in 1:4) {
    pdt <- pd[which(as.numeric(pd$fsCat)==(i-1)),]
    points(pval~d10.17,data=pdt,col=fsCols[i],type='b')
  }
  
}

fitARCMANmod <- function()
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  # model fitting - ARCMAN
  head(d)
  table(d$Species)
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  table(d$fsCat2) 
  
  fit5 <- glmmTMB(gCrown.18~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(gCrown.18~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(gCrown.18~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(gCrown.18~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (FALSE) 
  {
    fit5s <- brm(gCrown.18 ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    saveRDS(fit5s,'results/brm-ARCMAN-gc.RDS')
  } else {
    fit5s <- readRDS('results/brm-ARCMAN-gc.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F)
  P50r[1,1:6] <- c('ARCMAN','gCrown',p50)
  p50 <- plotPredictedValuesFit5(fit5s,T)
  P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitAMOCALmod <- function(yvar='gCrown')
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # model fitting - AMOCAL
  head(d)
  table(d$Species)
  table(d$SizeCat,d$fsCat)
  d$fsCat[which(d$fsCat==3)] <- 2
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  table(d$fsCat2) 
  
  if (yvar=='gCrown') d$yval <- d$gCrown.18 else d$yval <- d$DR.18
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yval),]
  dim(d)
  
  fit5 <- glmmTMB(yval~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yval~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yval~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yval~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (FALSE) 
  {
    fit5s <- brm(yval ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    if (yvar=='gCrown') saveRDS(fit5s,'results/brm-AMOCAL-gc.RDS') else saveRDS(fit5s,'results/brm-AMOCAL-dr.RDS')
  } else {
    if (yvar=='gCrown') fit5s <- readRDS('results/brm-AMOCAL-gc.RDS') else fit5s <- readRDS('results/brm-AMOCAL-dr.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  (p50 <- plotPredictedValuesFit5(fit5,F,yvar))
  P50r[1,1:6] <- c('AMOCAL',yvar,p50)
  (p50 <- plotPredictedValuesFit5(fit5s,T,yvar))
  P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitARBMENmod <- function(yvar='gCrown')
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # model fitting - ARBMEN
  head(d)
  table(d$Species)
  table(d$SizeCat,d$fsCat)
  d$SizeCat <- d$SSizeCat
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  table(d$SizeCat,d$fsCat)
  
  if (yvar=='gCrown') d$yval <- d$gCrown.18 else d$yval <- d$DR.18
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yval),]
  dim(d)
  
  fit5 <- glmmTMB(yval~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yval~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yval~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yval~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (FALSE) 
  {
    fit5s <- brm(yval ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    if (yvar=='gCrown') saveRDS(fit5s,'results/brm-ARBMEN-gc.RDS') else saveRDS(fit5s,'results/brm-ARBMEN-dr.RDS')
  } else {
    if (yvar=='gCrown') fit5s <- readRDS('results/brm-ARBMEN-gc.RDS') else fit5s <- readRDS('results/brm-ARBMEN-dr.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5_for_lrt3,F,yvar)
  P50r[1,1:6] <- c('ARBMEN',yvar,p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,yvar)
  P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitPSEMENmod <- function(logt=T)
{
  # model fitting - PSEMEN
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  head(d)
  table(d$Species)
  d$fsCat2 <- factor(d$fsCat)
  d$fsCat <- d$fsCat2
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  
  fit5 <- glmmTMB(gCrown.18~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(Live.18~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(Live.18~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(Live.18~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (FALSE) 
  {
    fit5s <- brm(gCrown.18 ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    saveRDS(fit5s,'results/brm-PSEMEN-gc.RDS')
  } else {
    fit5s <- readRDS('results/brm-PSEMEN-gc.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F,'gCrown',logt=logt)
  p50
  P50r[1,1:6] <- c('PSEMEN','gCrown',p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,'gCrown',logt=logt)
  p50
  P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitQUEAGRmod <- function(yvar='gCrown',logt=F,inclPlot=T)
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # model fitting - QUEAGR
  head(d)
  table(d$Species)
  table(d$SizeCat,d$fsCat)
  d$SizeCat <- d$TSizeCat
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  table(d$SizeCat,d$fsCat)
  
  if (yvar=='gCrown') d$yval <- d$gCrown.18 else if (yvar=='DR') d$yval <- d$DR.18 else d$yval <- d$Live.18
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yval),]
  dim(d)
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  
  fit5 <- glmmTMB(yval~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yval~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yval~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yval~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (FALSE) 
  {
    fit5s <- brm(yval ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   #(1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    if (yvar=='gCrown') saveRDS(fit5s,'results/brm-QUEAGR-gc.RDS') else if (yvar=='DR') saveRDS(fit5s,'results/brm-QUEAGR-dr.RDS') else saveRDS(fit5s,'results/brm-QUEAGR-lv.RDS')
  } else {
    if (yvar=='gCrown') fit5s <- readRDS('results/brm-QUEAGR-gc.RDS') else if (yvar=='DR') fit5s <- readRDS('results/brm-QUEAGR-dr.RDS') else readRDS('results/brm-QUEAGR-lv.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F,yvar)
  p50
  P50r[1,1:6] <- c('QUEAGR',yvar,p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,yvar)
  p50
  P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitUMBCALmod <- function(yvar='gCrown',logt=F)
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # model fitting - UMBCAL
  head(d)
  table(d$Species)
  table(d$SizeCat,d$fsCat)
  d$SizeCat <- d$TSizeCat
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  table(d$SizeCat,d$fsCat)
  
  if (yvar=='gCrown') d$yval <- d$gCrown.18 else if (yvar=='DR') d$yval <- d$DR.18 else d$yval <- d$Live.18
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yval),]
  dim(d)
  
  if (logt) d$d10.17 <- log10(d$d10.17)
  
  fit5 <- glmmTMB(yval~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yval~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yval~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yval~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (TRUE) 
  {
    fit5s <- brm(yval ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    if (yvar=='gCrown') saveRDS(fit5s,'results/brm-UMBCAL-gc.RDS') else if (yvar=='DR') saveRDS(fit5s,'results/brm-UMBCAL-dr.RDS') else saveRDS(fit5s,'results/brm-UMBCAL-lv.RDS')
  } else {
    if (yvar=='gCrown') fit5s <- readRDS('results/brm-UMBCAL-gc.RDS') else if (yvar=='DR') fit5s <- readRDS('results/brm-UMBCAL-dr.RDS') else readRDS('results/brm-UMBCAL-lv.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F,yvar)
  p50
  #P50r[1,1:6] <- c('UMBCAL',yvar,p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,yvar)
  p50
  #P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitQUEGARmod <- function(yvar='gCrown')
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # model fitting - QUEGAR
  head(d)
  table(d$Species)
  table(d$SizeCat,d$fsCat)
  d$SizeCat <- d$TSizeCat
  d$fsCat[which(d$fsCat==3)] <- 2
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  d$fsCat <- d$fsCat2
  table(d$SizeCat,d$fsCat2)
  
  if (yvar=='gCrown') d$yval <- d$gCrown.18 else if (yvar=='DR') d$yval <- d$DR.18 else d$yval <- d$Live.18
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yval),]
  table(d$SizeCat,d$fsCat)
  dim(d)
  
  fit5 <- glmmTMB(yval~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yval~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yval~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yval~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (TRUE) 
  {
    fit5s <- brm(yval ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    if (yvar=='gCrown') saveRDS(fit5s,'results/brm-QUEGAR-gc.RDS') else if (yvar=='DR') saveRDS(fit5s,'results/brm-QUEGAR-dr.RDS') else saveRDS(fit5s,'results/brm-QUEGAR-lv.RDS')
  } else {
    if (yvar=='gCrown') fit5s <- readRDS('results/brm-QUEGAR-gc.RDS') else if (yvar=='DR') fit5s <- readRDS('results/brm-QUEGAR-dr.RDS') else readRDS('results/brm-QUEGAR-lv.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F,yvar)
  p50
  P50r[1,1:6] <- c('QUEGAR',yvar,p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,yvar)
  p50
  P50r[1,7:10]<- p50
  P50r
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}

fitHETARBmod <- function(yvar='gCrown')
{
  P50r <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # model fitting - HETARB
  head(d)
  table(d$Species)
  table(d$SizeCat,d$fsCat)
  d$SizeCat <- d$SSizeCat
  d$fsCat2 <- factor(as.numeric(d$fsCat)-1)
  d$fsCat <- d$fsCat2
  table(d$SizeCat,d$fsCat)
  
  if (yvar=='gCrown') d$yval <- d$gCrown.18 else if (yvar=='DR') d$yval <- d$DR.18 else d$yval <- d$Live.18
  d <- d[complete.cases(d$fsCat2,d$d10.17,d$yval),]
  dim(d)
  
  fit5 <- glmmTMB(yval~d10.17 * fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  print(summary(fit5))
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(yval~d10.17 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt)
  print('fire severity effect')
  print(anova(fit5, fit5_for_lrt, "LRT"))
  
  fit5_for_lrt2 <- glmmTMB(yval~fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt2)
  print('size effect')
  print(anova(fit5, fit5_for_lrt2, "LRT"))
  
  fit5_for_lrt3 <- glmmTMB(yval~d10.17 + fsCat2 + (1|Plot) + (1|TreeNum), data=d, family='binomial')
  summary(fit5_for_lrt3)
  print('severity*size interaction effect')
  print(anova(fit5, fit5_for_lrt3, "LRT"))
  
  if (TRUE) 
  {
    fit5s <- brm(yval ~ 
                   s(d10.17, by=fsCat2, k=20)+
                   fsCat2+
                   (1|Plot)+
                   (1|TreeNum),
                 data=d,
                 family="Bernoulli",
                 chains = 2, cores = 2, seed=237, 
                 #backend="cmdstanr",
                 control=list(adapt_delta=0.99))
    conditional_effects(fit5s)
    if (yvar=='gCrown') saveRDS(fit5s,'results/brm-HETARB-gc.RDS') else if (yvar=='DR') saveRDS(fit5s,'results/brm-HETARB-dr.RDS') else saveRDS(fit5s,'results/brm-HETARB-lv.RDS')
  } else {
    if (yvar=='gCrown') fit5s <- readRDS('results/brm-HETARB-gc.RDS') else if (yvar=='DR') fit5s <- readRDS('results/brm-HETARB-dr.RDS') else readRDS('results/brm-HETARB-lv.RDS')
  }
  
  # second argument is F for glmmTMB models and T for brm models
  p50 <- plotPredictedValuesFit5(fit5,F,yvar)
  p50
  #P50r[1,1:6] <- c('HETARB',yvar,p50)
  p50 <- plotPredictedValuesFit5(fit5s,T,yvar)
  p50
  #P50r[1,7:10]<- p50
  
  return(P50r)
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
}
