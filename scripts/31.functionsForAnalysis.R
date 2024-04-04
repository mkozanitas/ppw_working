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
    barplot(table(temp$SizeCat,temp$fsCat,temp$fate3.18)[,,1]/tot,beside = T,col=tree.cols,main=paste(spname,'Mortality'),ylim=c(0,1))
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