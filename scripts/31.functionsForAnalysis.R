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
  
  TernaryPlot(alab='CON',blab='EHRO',clab='WHTO')
  TernaryText(pSGBA3,substr(row.names(pSGBA3),6,7))
  
  # classify plots
  pt <- data.frame(plot=row.names(pSGBA),vt=NA)
  pt$vt[which(pSGBA3[,1]>=50)] <- 'CON'
  pt$vt[which(pSGBA3[,2]>=50)] <- 'MH'
  pt$vt[which(pSGBA3[,3]>=50)] <- 'WO'
  #pt$vt[which(is.na(pt$vt) & pSGBA3[,1]<20)] <- 'MH-WO'
  #pt$vt[which(is.na(pt$vt) & pSGBA3[,3]<20)] <- 'MH-CON'
  pt$vt[which(is.na(pt$vt))] <- 'Mix3'
  pt
  write.csv(pt,'data/vegtypes.csv')
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

barplotNonSprouters <- function(print.to.pdf=T)
{
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  tdat <- list()
  
  if (print.to.pdf) pdf(paste('results/fates-non-sprouters.pdf',sep=''),6,6)

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
  if (print.to.pdf) dev.off()
  return(tdat)
}

barplotSprouters <- function(print.to.pdf=T)
{
  tree.cols <- c('grey90','grey60','grey30','black')
  shrub.cols <- c('grey90','black')
  tdat <- list()
  
  tsets <- list()
  tset.names <- c()
  tsets[[1]] <- c('QUEDOU','QUEGAR')
  tset.names[1] <- 'White Oaks'
  tsets[[2]] <- c('ARBMEN','QUEAGR','QUEKEL','UMBCAL')
  tset.names[2] <- 'EHRO'
  
  if (print.to.pdf) pdf(paste('results/fates-sprouters.pdf',sep=''),10,8)
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
  if (print.to.pdf) dev.off()
  return(tdat)
}