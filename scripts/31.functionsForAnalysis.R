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
  pt$vt[which(pSGBA3[,1]>70)] <- 'CON'
  pt$vt[which(pSGBA3[,2]>70)] <- 'MH'
  pt$vt[which(pSGBA3[,3]>70)] <- 'WO'
  pt$vt[which(is.na(pt$vt) & pSGBA3[,1]<20)] <- 'MH-WO'
  pt$vt[which(is.na(pt$vt) & pSGBA3[,3]<20)] <- 'MH-CON'
  pt$vt[which(is.na(pt$vt))] <- 'Mix3'
  write.csv(pt,'data/vegtypes.csv')
  return(pt)
}

calcFatesTableBySpecies <- function()
{
  fst12 <- data.frame(SpCode=rep(use.species,each=2),Type=rep(c('SA','TR'),length(use.species)),N17=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
  head(fst12)                  
  tail(fst12)
  
  i=56
  for (i in 1:nrow(fst12))
  {
    sp <- fst12$SpCode[i]
    ty <- fst12$Type[i]
    temp <- tAll[which(tAll$Species.17==sp & tAll$Type.17==ty),]
    
    fst12$N17[i] <- sum(temp$Live.17,na.rm=T)
    
    ## The next three lines are all equivalent - just using third one
    #fst12$N18.DN[i] <- length(which(temp$fate.18=='DN'))
    #fst12$N18.DN[i] <- length(which(temp$DN.18=='1'))
    fst12$N18.DN[i] <- sum(temp$DN.18,na.rm = T)
    
    fst12$N18.DR[i] <- sum(temp$DR.18,na.rm = T)
    fst12$N18.LN[i] <- sum(temp$LN.18,na.rm = T)
    fst12$N18.LR[i] <- sum(temp$LR.18,na.rm = T)
    miss <- which(temp$Live.13==1 & is.na(temp$DN.18)==1)
    if (length(miss)>0) for (j in 1:length(miss)) print(temp[miss[j],c('Plot.13','Num')])
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