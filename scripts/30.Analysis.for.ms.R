## This is final analysis script with code in order required to reproduce results in Results section of manuscript

rm(list=ls())
#require(MuMln)

source('scripts/31.functionsForAnalysis.R')
op.reset <- par(mfcol=c(1,1),mar=c(5,5,3,3))
#fsCols <- c('blue','brown','orange','red')
fsCols <- c('#00FF01','#FFB000','#FE6100','#DC267F')
fsCols3 <- c("0.U" = '#00FF01',
        "12.LM" = '#FE8800',
        "3.H" = '#DC267F')
vtShapes <- data.frame(vt=c('CON','MH','Mix3','WO'),shp=c(18,15,17,16))

fs <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodFireSeverity/master/data/FSextract/vegplots-54-20m-FS.csv")
head(fs)
tail(fs)
dim(fs)

hfs <- read.csv('input_data/plot_info/hectares-18-20m-FS.csv')
hfs$Plot.Orig <- hfs$Plot
substr(hfs$Plot,4,5) <- '13'
names(hfs)[which(names(hfs)=='Quad')] <- 'Subplot'
head(hfs)

fs <- loadFireSeverity(fs,'Tubbs.MTBS.RDNBR.30',verbose=T)
head(fs)
head(fs[,c('Plot','fsvar','fsCat')])
table(fs$fsCat)

hfs <- loadFireSeverity(hfs,'Tubbs.MTBS.RDNBR.30',verbose=F)
head(hfs)
table(hfs$fsCat)

if (FALSE) {
  cols <- c('purple','lightblue','darkblue','red')
  plot(fs$Tubbs.MTBS.RDNBR.30,fs$Kincade.MTBS.RDNBR.30,pch=19,cex=2,xlab='Tubbs Fire',ylab='Kincade Fire',main='MTBS-dNBR',col=cols[fs$fsLevel+1])
  abline(v=c(135,430),lty=2)
  abline(h=c(135,430),lty=2)
}

# Read in species codes - in this script, called 'Species'
spAtt <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/SpeciesData/all-spp-attributes.csv")

dim(spAtt)
head(spAtt)
tail(spAtt)
table(spAtt$FuncGroup)

# input merged dataframe - read in with hectare data
tAll <- read.csv('data/tAllh.csv',as.is=T,row.names=1)
head(tAll)
dim(tAll)
names(tAll)
table(tAll$Survey)

# convert Species# convert QUEBEGA to QUEBER, QUEDOGA to QUEDOU, QUEdec to QUEDOU, and remove unknowns
table(tAll$Species)
tAll <- convertHybrids()
table(tAll$Species)
sum(table(tAll$Species))

# add functional groups to tAll
# Are all species in spAtt?
(tAllspec <- sort(unique(tAll$Species)))
tAllspec[which(!tAllspec %in% spAtt$OrigSpecies)]
length(which(tAll$Species=='QUEDEC'))

tAll$FuncGroup <- spAtt$FuncGroup[match(tAll$Species,spAtt$OrigSpecies)]
table(tAll$FuncGroup)

# remove Type.13 = 'TS'; a few show up in 2018 - need to check
dim(tAll)
tAll <- tAll[-which(tAll$Type.17=='TS'),]
dim(tAll)

# assign plot and hectare rows
pRows <- which(tAll$Survey=='Plot')
hRows <- which(tAll$Survey=='Hect')

# summary statistics for fire severity by plots
table(fs$fsCat)
table(hfs$fsCat,hfs$Plot)
table(hfs$fsCat)

# add fire severity to tAll
table(tAll$Plot,useNA='always')
table(tAll$Plot.Orig,useNA='always')
fs$Plot.Orig <- fs$Plot
substr(fs$Plot,4,5) <- '13'
tAll$fsCat <- NA
tAll$fsCat[pRows] <- fs$fsCat[match(tAll$Plot[pRows],fs$Plot)]
table(tAll$fsCat[pRows],useNA='always')

tAll$PSp <- NA
tAll$PSp[hRows] <- paste(tAll$Plot[hRows],tAll$Subplot.17[which(tAll$Survey=='Hect')])
hfs$PSp <- paste(hfs$Plot,hfs$Subplot)
length(which(!tAll$PSp[hRows] %in% hfs$PSp))
tAll$fsCat[hRows] <- hfs$fsCat[match(tAll$PSp[hRows],hfs$PSp)]
table(tAll$fsCat[hRows],useNA='always')
table(tAll$fsCat,useNA='always')

# archive tAll in tAll.arch and then reduce tAll to allow for easier examination during analysis
tAll.archive <- tAll
table(tAll.archive$Species)

# use this to restore and recreate tAll - create dataframes with and without hectares
tAll <- tAll.archive
table(tAll$ExcStem)

tAll$BABH.17 <- pi*(tAll$DBH_cm.17/2^2)
tAll$BAd10.17 <- pi*(tAll$d10.17/2^2)

# subset to live prefire - removes individuals newly added in 2019 and 2020
tAll <- tAll[which(tAll$Live.17==1),c('Num','TreeNum','Species','FuncGroup','Plot','Subplot.17','fsCat','Year.13','Point.13','Year.17','Type.17','Live.17','DBH_cm.17','d10.17','BABH.17','BAd10.17','SA.Height_cm.17','Year.18','Type.18','Live.18','fate3.18','DN.18','DR.18','LN.18','LR.18','gCrown.18','DBH_cm.18','d10.18','Basal.Resprout.Count.18','Basal.Resprout.Height_cm.18','Year.19','Live.19','fate3.19','Type.19','DBH_cm.19','d10.19','Basal.Resprout.Count.19','Basal.Resprout.Height_cm.19','Survey')]

# Remove Baccharis and QUEDEC
tAll <- tAll[-which(tAll$Species %in% c('BACPIL','QUEDEC')),]
dim(tAll)
write.csv(tAll,"data/tAll.csv")

table(tAll$fate3.18,tAll$Survey,useNA='a')
tAll <- tAll[-which(is.na(tAll$fate3.18)),]
table(tAll$fate3.18,tAll$Survey,useNA='a')

# Size distribution in four levels
tAll$Type.17.4 <- NA
tAll$Type.17.4[which(tAll$Live.17==1 & tAll$Type.17=='SA')] <- 'SA'
tAll$Type.17.4[which(tAll$Live.17==1 & tAll$Type.17=='TR' & tAll$DBH_cm.17<= 10)] <- 'ST'
tAll$Type.17.4[which(tAll$Live.17==1 & tAll$Type.17=='TR' & tAll$DBH_cm.17> 10 & tAll$DBH_cm.17<=20)] <- 'MT'
tAll$Type.17.4[which(tAll$Live.17==1 & tAll$Type.17=='TR' & tAll$DBH_cm.17> 20)] <- 'LT'

table(tAll$Type.17.4,useNA='a')
summary(tAll$DBH_cm.17[which(is.na(tAll$Type.17.4))])
summary(tAll$d10.17[which(is.na(tAll$Type.17.4))])
table(tAll$Type.17[which(is.na(tAll$Type.17.4))])

# now check d10 to try and assign a few remaining individuals
tAll$Type.17.4[which(is.na(tAll$Type.17.4) & tAll$Type.17=='TR' & tAll$d10.17 <= 12.83)] <- 'ST'
tAll$Type.17.4[which(is.na(tAll$Type.17.4) & tAll$Type.17=='TR' & tAll$d10.17 > 12.83 & tAll$d10.17<=24.59)] <- 'MT'
tAll$Type.17.4[which(is.na(tAll$Type.17.4) & tAll$Type.17=='TR' & tAll$d10.17> 24.59)] <- 'LT'
tAll$Type.17.4[which(tAll$Type.17.4=='LT' & tAll$Survey=='Hect')] <- 'HLT'
table(tAll$Type.17.4,useNA='a')
table(tAll$fsCat,tAll$Type.17.4,useNA='a')

## Drop 37 individuals that are not classified by size and 37 individuals (coincidence it's the same?) without FS
tAll <- tAll[-which(is.na(tAll$Type.17.4)),]
tAll <- tAll[-which(is.na(tAll$fsCat)),]

# Table 2
(T2 <- table(tAll$Type.17.4,tAll$fsCat,useNA='a'))
sum(T2)

# Table 1
T1 <- table(tAll$Species)
prows <- which(tAll$Survey=='Plot')
hrows <- which(tAll$Survey=='Hect')
T1N.p <- table(tAll$Species[prows])
T1N.p <- T1N.p[match(names(T1),names(T1N.p))]
T1N.h <- table(tAll$Species[hrows])
T1N.h <- T1N.h[match(names(T1),names(T1N.h))]
T1babh.plots <- tapply(tAll$BABH.17[prows],tAll$Species[prows],sum,na.rm=T)
T1babh.plots <- T1babh.plots[match(names(T1),names(T1babh.plots))]
T1babh.hect <- tapply(tAll$BABH.17[hrows],tAll$Species[hrows],sum,na.rm=T)
T1babh.hect <- T1babh.hect[match(names(T1),names(T1babh.hect))]
T1babh <- tapply(tAll$BABH.17,tAll$Species,sum,na.rm=T)
T1bad10 <- tapply(tAll$BAd10.17,tAll$Species,sum,na.rm=T)
T1 <- cbind(T1,T1N.p,T1N.h,T1babh.plots,T1babh.hect,T1bad10)
T1
T1 <- T1[order(T1[,1],decreasing = T),]
write.csv(T1,'results/table1.csv')
summary(tAll$BABH.17)
summary(tAll$BAd10.17)

# THIS IS FINAL ANALYSIS SET
# Introduction - lists number of stems sampled. Include 2018 plots as this is just a general statement about size of the network
dim(tAll)
table(tAll$Survey,useNA='a')

# at this point, can subset to Plot or Hect surveys
pRows <- which(tAll$Survey=='Plot')
hRows <- which(tAll$Survey=='Hect')

# add Subplot to Plot data
tAll$Subplot.17[pRows] <- 'C3'
table(tAll$Subplot.17,tAll$Survey)

# mortality in unburned, as control on rates of mortality we may have missed before fire
unb.mort <- table(tAll$Type.17[which(tAll$fsCat==1)],tAll$fate3.18[which(tAll$fsCat==1)])
sum(unb.mort[,1])/sum(unb.mort)

# RESULTS
table(fs$fsCat)

# First result - ternary plot ordination of plots to describe veg types, intersected with fire severity
pt <- drawTernaryPlots(d=tAll[pRows,])

head(pt)
tail(pt)
if (all(pt$plot==fs$Plot)) pt$fsCat <- fs$fsCat else print('error')
head(pt)
tail(pt)
table(pt$vt)
write.csv(pt,'results/plot-types-fs.csv')
table(pt$fsLevel)

(vtfs <- table(pt[,c('vt','fsCat')]))
write.csv(vtfs,'results/veg-type-fire-sev-table.csv')

# Second result - overall numbers
dim(tAll)

# total number of large trees
length(which(tAll$Type.17.4 %in% c('LT','HLT')))

# count by tree size and whether in plots of hectares
table(tAll$Type.17.4,tAll$Survey)

# TABLE B (2)
table(fs$fsCat)
table(hfs$fsCat)
table(tAll$Type.17.4,tAll$fsCat)

# summary of fates - All SA and TR, PLOTS
(ftab <- table(tAll$fate3.18[pRows]))
(allAb18 <- sum(table(tAll$fate3.18[pRows])))
ftab/allAb18 # High level of 2018 fates

# summary of fates - Large Trees
ltRows <- which(tAll$Type.17.4 %in% c('LT','HLT'))
(ftab <- table(tAll$fate3.18[ltRows]))
(allAb18 <- sum(table(tAll$fate3.18[ltRows])))
ftab/allAb18 # High level of 2018 fates

# MISCELLANEOUS STATS NOT SHOWN IN MS CURRENTLY
if (FALSE) {
  # Number of live species prefire in plots
  table(tAll$Species)
  length(unique(tAll$Species))
  length(unique(tAll$Species[which(tAll$Survey=='Plot')]))
  length(unique(tAll$Species[which(tAll$Survey=='Hect')]))
  table(tAll$Species,tAll$Survey)
  
  table(tAll$Type.17,tAll$Survey)
  sum(table(tAll$Type.17,tAll$Survey))
  table(tAll$Type.17[ltRows],tAll$Survey[ltRows])
  
  #Species abundance in PLOTS (hect removed to avoid mixing size class data)
  (allAb <- sum(table(tAll$Type.17[pRows])))
  (spAbund <- sort(table(tAll$Species[pRows]),decreasing = T))
  (spAbund.all <- table(tAll$Species))
  a2s <- match(spAtt$OrigSpecies,names(spAbund.all))
  spAtt$tot.abund <- spAbund.all[a2s]
  write.csv(spAtt,'results/all-spp-attributes.csv')
  
  # PSEMEN as proportion of all non-sprouters
  spAbund[which(names(spAbund)=='PSEMEN')]
  spAbund[which(names(spAbund)=='PSEMEN')]/sum(spAbund[which(names(spAbund) %in% spAtt$Species[which(spAtt$Resprout=='N')])])
  
  #### Summary stats
  (use.species <- spAtt$Species)
  
  (Nresprouters <- sum(spAbund[which(names(spAbund) %in% spAtt$Species[which(spAtt$Resprout=='Y')])]))
  
  (common.species <- spAtt$Species[which(spAtt$Common=='Yes')])
  (sumAb <- sum(spAbund[which(names(spAbund) %in% common.species)]))
  sumAb/allAb
  
  ## Abundance of common species for PLOTS + HECTS
  (spAbund.ph <- sort(table(tAll$Species),decreasing = T))
  (spAbund.ph.common <- sort(table(tAll$Species[which(tAll$Species %in% common.species)]),decreasing = T))
  sum(spAbund.ph.common)/sum(spAbund.ph)
  
  
  # MELINA: QUESTION
  if (FALSE) {
    # how many 'new' saplings 
    newSap <- which(is.na(tAll$Year.13) & tAll$Year.18==2018 & tAll$Type.18=='SA')
    table(tAll$Plot[newSap])
    length(newSap)
        
    hist(tAll$SA.Height_cm.17[newSap])
    newTr <- which(is.na(tAll$Year.13) & tAll$Year.18==2018 & tAll$Type.18=='TR')
    
    # how many 'new' trees
    #tAllm$Num[newSap]
    #table(tAllm$Plot[newSap])
    #table(tAllm$fate.18[newSap],tAllm$fsCat[newSap])
    
    # MELINA: Remind me about these - they were either tagged but data failed to record in 2013, or newly recruited from 2013 to 2017 and inferred in '2017' data. They are detected here by NA in Year.17 which I think means 
  }
  
  # summary of species types
  table(spAtt$Shrub.Tree)
}

# create fst dataframe - FateSummaryTable for time 1 -> 2 (2017 and 2018) - PLOTS
use.species <- sort(unique(tAll$Species))
fst12 <- calcFatesTableBySpecies(use.species,survey=c('Plot','Hect'))
head(fst12)
tail(fst12)
write.csv(fst12,'results/fates.x.species.csv')

table(tAll$Type.17.4,tAll$fate3.18)
table(tAll$Type.17.4)
table(tAll$fate3.18)

#NOT USED IN MS
if (FALSE) {
  # reduce fst12 table to common species plus other shrubs and trees each tallied by SA and TR
  fst12c <- reduce_fst12()
  fst12c
  for (i in 3:7) fst12c[,i] <- as.numeric(fst12c[,i])
  
  # column sums of fst12c, and as percentage of total
  (allSum <- apply(as.matrix(fst12c[,-c(1:2)]),2,sum))
  allSum/allSum[1]
  
  rs <- which(fst12c$Type=='SA')
  (SASum <- apply(as.matrix(fst12c[rs,-c(1:2)]),2,sum))
  rs <- which(fst12c$Type=='TR')
  (TRSum <- apply(as.matrix(fst12c[rs,-c(1:2)]),2,sum))
  
  ## summary of fates for all, saplings alone, trees alone
  allSum/allSum[1]
  SASum/SASum[1]
  TRSum/TRSum[1]
  
  write.csv(fst12c,'results/summary-table.csv')
  
  # preliminary poster results, not used in the end
  if (FALSE) {
    poster.species <- c('PSEMEN','QUEAGR','UMBCAL','QUEGAR','HETARB')
    fst12p <- fst12c[which(fst12c$SpCode %in% poster.species),]
    fst12p
    
    fst12p[,4:7]/fst12p[,3]
    
    (allSum.p <- apply(as.matrix(fst12p[,-c(1:2)]),2,sum))
    allSum.p/allSum.p[1]
    rs <- which(fst12p$Type=='SA')
    (SASum.p <- apply(as.matrix(fst12p[rs,-c(1:2)]),2,sum))
    rs <- which(fst12p$Type=='TR')
    (TRSum.p <- apply(as.matrix(fst12p[rs,-c(1:2)]),2,sum))
    
    allSum.p/allSum.p[1]
    SASum.p/SASum.p[1]
    TRSum.p/TRSum.p[1]
    
    allSum.p/allSum
  }
}

# FOR NOW NOT UPDATING BARPLOTS FOR MS
if (FALSE) {
  # Barplots for sprouters and non-sprouters, all species
  fbp.list <- barplotFates(tAll,fs='low-medium')
  
  rspd <- tAll[which(tAll$Species %in% spAtt$Species[which(spAtt$Resprout=='Y')]),]
  fbp.rsp.list <- barplotFates(rspd,fs='low-medium')
  table(rspd$Type.17)
  
  nspd <- tAll[which(tAll$Species %in% spAtt$Species[which(spAtt$Resprout=='N')]),]
  fbp.nsp.list <- barplotFates(nspd,fs='low-medium')
  fbp.nsp.list
  table(nspd$Type.17)
  
  # non sprouters as proportion of all individuals
  table(nspd$Type.17)/table(tAll$Type.17)
  sum(table(nspd$Type.17))/sum(table(tAll$Type.17))
  
  #### makeffsp2#### make bar plotsum#### make bar plots for different subgroups and individual species
  # tAll has all the living plants from 2017, with their 2018 fates
  # doesn't have new plants added in 2019
  uh <- TRUE #whether to use hectares for model fitting
  if (uh) dat <- tAll else dat = tAll[pRows,]
  tAllm <- prepareForBarPlots(d=dat)
  dim(tAll)
  dim(tAllm)
  table(tAllm$TSizeCat)
  sum(table(tAllm$TSizeCat))
  sort(table(tAllm$Species),decreasing = T)
  
  # results holder for P50 results - critical size to reach 50% green crown or topkill.resprout, based on either maxlikelihood or baseyian models
  #P50res <- data.frame(Species=NA,yvar=NA,GCml.0=NA,GCml.1=0,GCml.2=0,GCml.3=0,GCby.0=NA,GCby.1=0,GCby.2=0,GCby.3=0)
  
  # barplot for PSEMEN, first to see here, and then to pdf
  tdat <- barplotNonSprouters(d=tAllm)
  length(tdat)
  d <- tdat[[1]]
  table(d$Species)
}

## Next three sections cover PSEMEN, then all the resprouters, and then the resprouting functional groups. These are focused on building the models and the figures showing fate3 as a function of size and fire severity, by species or FG. 

#In order to have the same range of x-axes on all figures, this next line identifying the range of size values. Then for each figure it's adjusted to only calculate and show the curves for the span of sizes covered by that species or group
drange <- range(tAll$d10.17,na.rm=T)
log10(drange)

spSel <- 'PSEMEN'
spName <- spSel
fs=c('low-medium') #'all','low-medium'
#fs=c('drop-high','low-medium') #AMOCAL, QUEGAR
#fs=c('low-medium','drop-high','drop-unburned') #FRACAL
iter=10000
logt=T
uh <- TRUE

if (uh) d <- tAll[which(tAll$Species == spSel),] else d <- tAll[which(tAll$Species == spSel & tAll$Survey=='Plot'),]
range(log10(d$d10.17),na.rm=T)
dim(d)
names(d)
table(d$Survey)

table(d$Species)

table(as.numeric(d$fsCat))
fac.fsCat.levels <- c('0.U','1.L','2.M','3.H')
d$fac.fsCat <- fac.fsCat.levels[as.numeric(d$fsCat)]
table(d$fac.fsCat,d$Survey)

fslevels <- c()
if ('all' %in% fs) {
  fslevels <- c(fslevels,'fs.all')
}
if ('drop-high' %in% fs)
{
  d$fac.fsCat[which(d$fac.fsCat=='3.H')] <- '2.M'
  fslevels <- c(fslevels,'drop.high')
}
if ('drop-unburned' %in% fs)
{
  d$fac.fsCat[which(d$fac.fsCat=='0.U')] <- '1.L'
  fslevels <- c(fslevels,'drop-unburned')
}
if ('low-medium' %in% fs)
{
  d$fac.fsCat[which(d$fac.fsCat %in% c('1.L','2.M'))] <- '12.LM'
  fslevels <- c(fslevels,'comb-low-med')
}
print(fslevels)
table(d$fac.fsCat,d$Survey)
d$Live.18[which(d$fate3.18=='GC')] <- 1
d$Live.18[which(d$fate3.18=='DN')] <- 0
table(d$Live.18,d$fate3.18,useNA='a')

dd <- d[complete.cases(d$fac.fsCat,d$d10.17,d$Live.18,d$Plot,d$TreeNum),]
table(dd$Survey)
if (logt) dd$d10.17 <- log10(dd$d10.17)
table(d$fate3.18,d$fac.fsCat)
dim(dd)
names(dd)

saveRDS(dd,paste(local.dir,'/brm.',spName,'.dd.rds',sep=''))

table(dd$Plot,dd$fac.fsCat)
spSel

tdat <- barplotOneNonSprouter(dd,spSel,T,T,F)

# model fitting, if needed 
if (FALSE) {
  fit5brm <- brm(Live.18 ~ s(d10.17, k=3, by=fac.fsCat) + fac.fsCat + (1|Plot), data=dd,
               family="bernoulli", 
               chains = 2,
               cores = 2, 
               seed=726, 
               iter=iter,
               #backend="cmdstanr",
               refresh=100,
               control=list(adapt_delta=0.99))
  # fit5brm.noplot <- brm(Live.18 ~ s(d10.17, k=3, by=fac.fsCat) + fac.fsCat, data=dd,
  #                family="bernoulli", 
  #                chains = 2,
  #                cores = 2, 
  #                seed=726, 
  #                iter=iter,
  #                #backend="cmdstanr",
  #                refresh=100,
  #                control=list(adapt_delta=0.99))
  beep()
saveRDS(warnings(),paste(local.dir,'/brm.',spName,'.',uh,'.',iter,'.BERN.Splk3.Live18.WARNINGS.rds',sep=''))
print(summary(fit5brm))
saveRDS(fit5brm,paste(local.dir,'/brm.',spName,'.',uh,'.',iter,'.BERN.Splk3.Live18.rds',sep=''))

multifit <- readRDS(paste(local.dir,'/brm.',spName,'.',uh,'.',iter,'.BERN.Splk3.Live18.rds',sep=''))
visualizeBernfitBayes(mf=multifit,sp='',print.to.pdf=T,xlims=c(-1,2.3),spName=spName)

multifit <- readRDS(paste(local.dir,'/brm.',spName,'.noplot.',uh,'.',iter,'.BERN.Splk3.Live18.rds',sep=''))
visualizeBernfitBayes(mf=multifit,sp='',print.to.pdf=T,xlims=c(-1,2.3),spName=spName)

}

# refresh script 31 if needed
source('scripts/31.functionsForAnalysis.R')

spList <- sort(table(tAll$Species),decreasing=T)
spList <- spList[-which(spList<400)]
(spList <- spList[-which(names(spList)=='PSEMEN')])

## Set species, fire severity option, and log-size option

## CODE FOR INDIVIDUAL SPECIES - RUN INTERACTIVELY
i=7

# BAR PLOT FUNCTION NOT WORKING FOR SHRUBS
for (i in 5:length(spList)) 
{
  spSel <- names(spList)[i]
  spName <- spSel
  fs=c('low-medium') #'all','low-medium'
  if (spSel %in% c('AMOCAL','QUEGAR','QUEKEL','QUEDOU')) fs=c('drop-high','low-medium')
  if (spSel %in% c('FRACAL')) fs=c('low-medium','drop-high','drop-unburned')
  #if (spSel %in% c('QUEGAR','QUEDOU','QUEKEL')) fs=c('low-medium','drop-high','drop-unburned') # analyze for burned only, with no FS term
  logt=T
  iter=50000
  uh <- TRUE
  
  
  #tdat <- barplotSprouterSpecies(spSel,skip.op=T)
  #dim(tdat)
  #tdat
  
  if (uh) d <- tAll[which(tAll$Species == spSel),] else d <- tAll[which(tAll$Species == spSel & tAll$Survey=='Plot'),]
  range(log10(d$d10.17),na.rm=T)
  dim(d)
  
  table(d$Species)
  
  table(as.numeric(d$fsCat))
  fac.fsCat.levels <- c('0.U','1.L','2.M','3.H')
  d$fac.fsCat <- fac.fsCat.levels[as.numeric(d$fsCat)]
  table(d$fac.fsCat)
  
  fslevels <- c()
  if ('all' %in% fs) {
    fslevels <- c(fslevels,'fs.all')
  }
  if ('drop-high' %in% fs)
  {
    d$fac.fsCat[which(d$fac.fsCat=='3.H')] <- '2.M'
    fslevels <- c(fslevels,'drop.high')
  }
  if ('drop-unburned' %in% fs)
  {
    d <- d[-which(d$fac.fsCat=='0.U'),]
    #d$fac.fsCat[which(d$fac.fsCat=='0.U')] <- '1.L'
    fslevels <- c(fslevels,'drop-unburned')
  }
  if ('low-medium' %in% fs)
  {
    d$fac.fsCat[which(d$fac.fsCat %in% c('1.L','2.M'))] <- '12.LM'
    fslevels <- c(fslevels,'comb-low-med')
  }
  print(fslevels)
  table(d$fac.fsCat)
  
  # fit multinomial first
  dd <- d[complete.cases(d$fac.fsCat,d$d10.17,d$fate3.18),]
  if (logt) dd$d10.17 <- log10(dd$d10.17)
  dim(dd)
  table(dd$fate3.18,dd$fac.fsCat)
  table(dd$fate3.18,dd$fac.fsCat,dd$Type.17.4)
  
  saveRDS(dd,paste(local.dir,'/brm.',spName,'.dd.rds',sep=''))
  tdat <- barplotSprouterSpecies(dd,spSel,print.to.pdf=T,T,F)
  
  # Run interactively
  # # this function runs k=3 spline, with more iterations
  # fitFatesMultinomial2.brm(d,spName,fs,logt,iter=2000)
  visualizeMultifitBayes_redraw(spName,sp='MN.Splk3')
  
}

## CODE FOR FUNCTIONAL GROUPS - RUN INTERACTIVELY
# now functional groups
# Functional Groups and sample sizes
#EHRO   NS.Con NS.Shrub  R.Shrub     WHTO 
#2850     1208      287     1845      477 

table(spAtt$FuncGroup,useNA='always')
FGroupList <- c('EHRO','WHTO','R.Shrub')
i=1
for (i in 1:length(FGroupList))
{
  FSel <- FGroupList[i]
  spName <- FSel
  (spSel <- spAtt$OrigSpecies[which(spAtt$FuncGroup==FSel)])
  if (FSel=='R.Shrub') spSel <- spSel[-which(spSel=='BACPIL')]
  sRow <- which(tAll$Species %in% spSel)
  table(tAll$Species[sRow])
  sum(table(tAll$Species[sRow]))
  
  fs=c('low-medium') #'all','low-medium'
  if (spName %in% c('WHTO','R.Shrub')) fs=c('drop-high','low-medium')
  #if (spName %in% c('WHTO','R.Shrub')) fs=c('low-medium','drop-high','drop-unburned') # analyze for burned only, with no FS term
  logt=T
  
  #tdat <- barplotSprouterSpecies(spSel,ss.name=spName,skip.op=T)
  #dim(tdat)
  
  if (uh) d <- tAll[which(tAll$FuncGroup %in% FSel),] else d <- tAll[which(tAll$Species == spSel & tAll$Survey=='Plot'),]
  range(log10(d$d10.17),na.rm=T)
  dim(d)
  
  table(d$Species)
  
  table(as.numeric(d$fsCat))
  fac.fsCat.levels <- c('0.U','1.L','2.M','3.H')
  d$fac.fsCat <- fac.fsCat.levels[as.numeric(d$fsCat)]
  table(d$fac.fsCat)
  
  fslevels <- c()
  if ('all' %in% fs) {
    fslevels <- c(fslevels,'fs.all')
  }
  if ('drop-high' %in% fs)
  {
    d$fac.fsCat[which(d$fac.fsCat=='3.H')] <- '2.M'
    fslevels <- c(fslevels,'drop.high')
  }
  if ('drop-unburned' %in% fs)
  {
    d <- d[-which(d$fac.fsCat=='0.U'),]
    #d$fac.fsCat[which(d$fac.fsCat=='0.U')] <- '1.L'
    fslevels <- c(fslevels,'drop-unburned')
  }
  if ('low-medium' %in% fs)
  {
    d$fac.fsCat[which(d$fac.fsCat %in% c('1.L','2.M'))] <- '12.LM'
    fslevels <- c(fslevels,'comb-low-med')
  }
  print(fslevels)
  table(d$fac.fsCat)
  
  # fit multinomial first
  dd <- d[complete.cases(d$fac.fsCat,d$d10.17,d$fate3.18),]
  if (logt) dd$d10.17 <- log10(dd$d10.17)
  dim(dd)
  table(dd$fate3.18,dd$fac.fsCat)
  
  saveRDS(dd,paste(local.dir,'/brm.',spName,'.dd.rds',sep=''))
  
  tdat <- barplotSprouterSpecies(dd,spSel,T,T,F,ss.name=FSel)
  
  # Run interactively
  # # this function runs k=3 spline, with more iterations
  # fitFatesMultinomial2.brm(d,spName,fs,logt,iter=2000)
  visualizeMultifitBayes_redraw(spName,sp='MN.Splk3')
  
}

if (FALSE) {
  # next four line pairs run 3 different spline models and then quadratic - I tested all of these, and settled on spline k=3, which is implemented above with more iterations to help with convergence.
  k=3 #3, 6, or 20 only
  fitFatesMultinomial.brm(d,spName,fs,logt,m.choice='spline',splk=k)
  
  k=6 #3, 6, or 20 only
  fitFatesMultinomial.brm(d,spName,fs,logt,m.choice='spline',splk=k)
  
  k=20 #3, 6, or 20 only
  fitFatesMultinomial.brm(d,spName,fs,logt,m.choice='spline',splk=k)
  
  fitFatesMultinomial.brm(d,spName,fs,logt,m.choice='quad')
}

## what about delayed mortality
# function only works for Plot data for now
delayedMortality(survey='Plot')

## resprout height growth
survey <- 'Plot'
# BASAL RESPROUTS script, in brackets
{
  tAll$Basal.Resprout.Height_cm.18 <- as.numeric(tAll$Basal.Resprout.Height_cm.18)
  tAll$Basal.Resprout.Height_cm.19 <- as.numeric(tAll$Basal.Resprout.Height_cm.19)
  tAll$Basal.Resprout.Count.18[which(tAll$Basal.Resprout.Count.18==0)] <- NA
  tAll$Basal.Resprout.Height_cm.18[which(tAll$Basal.Resprout.Height_cm.18==0)] <- NA
  tAll$Basal.Resprout.Count.19[which(tAll$Basal.Resprout.Count.19==0)] <- NA
  tAll$Basal.Resprout.Height_cm.19[which(tAll$Basal.Resprout.Height_cm.19==0)] <- NA
  
  # fix two outliers - remove this after data file is fixed
  tAll$Basal.Resprout.Height_cm.18[which.max(tAll$Basal.Resprout.Height_cm.18)] <- NA
  tAll$Basal.Resprout.Count.18[which.max(tAll$Basal.Resprout.Count.18)] <- NA
  
  # PSEMEN with basal resprouts 
  xx <- which(tAll$Basal.Resprout.Height_cm.18>=0 & tAll$Species=='PSEMEN')
  tAll[xx,]
  tAll$Basal.Resprout.Count.18[xx] <- NA
  tAll$Basal.Resprout.Height_cm.18[xx] <- NA
  xx <- which(tAll$Basal.Resprout.Height_cm.19>=0 & tAll$Species=='PSEMEN')
  tAll[xx,]
  tAll$Basal.Resprout.Count.19[xx] <- NA
  tAll$Basal.Resprout.Height_cm.19[xx] <- NA
  
  names(tAll)
  message('average resprout count and height - all species')
  print(mean(tAll$Basal.Resprout.Count.18,na.rm=T))
  print(mean(tAll$Basal.Resprout.Height_cm.18,na.rm=T))
  message('max resprout count and height')
  print(max(tAll$Basal.Resprout.Count.18,na.rm=T))
  print(max(tAll$Basal.Resprout.Height_cm.18,na.rm=T))
  
  cs <- spAtt$Species[which(spAtt$Shrub.Tree %in% c('S','T') & spAtt$Common=='Yes' & spAtt$Resprout=='Y')]
  tAllc <- tAll[which(tAll$Species %in% cs & (tAll$DR.18==1 | tAll$LR.18==1) & tAll$Survey %in% survey),]
  dim(tAllc)
  head(tAllc)
  
  tAllc$lBasal.Resprout.Count.18 <- log10(tAllc$Basal.Resprout.Count.18+1)
  hist(tAllc$lBasal.Resprout.Count.18)
  hist(log10(tAllc$Basal.Resprout.Height_cm.18))
  plot(tAllc$lBasal.Resprout.Count.18,log10(tAllc$Basal.Resprout.Count.19+1))
  abline(0,1)
  
  fit <- lm(Basal.Resprout.Height_cm.18~Species*fsCat,data=tAllc)
  anova(fit)
  summary(fit)
  hist(fit$residuals)
  par(mfrow=c(2,1),mar=c(5,5,3,1))
  boxplot(Basal.Resprout.Height_cm.18~Species,data=tAllc,ylab='Resprout height (cm), 2018')
  boxplot(Basal.Resprout.Height_cm.18~fsCat,data=tAllc,xlab='Fire Severity',ylab='Resprout height (cm), 2018')
  
  
  length(which(tAllc$Basal.Resprout.Height_cm.18>0 & tAllc$fsCat==1))
  
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
  par(mfrow=c(2,1),mar=c(5,5,3,1))
  boxplot(RespGrowth~Species,data=tAllc,ylab='Resprout growth, 2018-19 (cm)')
  abline(h=0,lty=2)
  boxplot(RespGrowth~fsCat,data=tAllc,xlab='Fire Severity',ylab='Resprout growth, 2018-19 (cm)')
  abline(h=0,lty=2)
  
  par(op.reset)
  plot(Basal.Resprout.Height_cm.19~Basal.Resprout.Height_cm.18,data=tAllc,xlab='Basal Resprout Height (cm), 2018',ylab='Basal Resprout Height (cm), 2019',asp=1,xlim=c(0,500),ylim=c(0,500),col=fsCols[as.numeric(tAllc$fsCat)],pch=19,cex=1.5,cex.lab=1.5)
  fit <- lm(Basal.Resprout.Height_cm.19~Basal.Resprout.Height_cm.18,data=tAllc)
  abline(fit)
  anova(fit)
  abline(0,1,lty=2)
  
  relGrowth <- tAllc$Basal.Resprout.Height_cm.19/tAllc$Basal.Resprout.Height_cm.18
  summary(relGrowth)
}

# 2019 recruits
length(which(tAll$Live.18==1))
length(which(tAll$Live.18==1 & tAll$Type.18=='SA'))


# NEW RECRUITS script, in brackets
# FIX THIS
{
  # restore from archive to get new individuals
  t19 <- tAll.archive
  t19 <- t19[-which(t19$Species %in% c('BACPIL','QUEDEC')),]
  t19 <- t19[which(is.na(t19$fate3.18) & !is.na(t19$fate3.19)),]
  t19 <- t19[-which(t19$Type.19=='TS'),]
  dim(t19)
  
  print('number of new recruits, by type, species, fire severity')
  table(t19$Type.19,useNA='a')
  table(t19$Species,t19$Type.19,useNA='a')
  
  t19tr <- which(t19$Type.19=='TR')
  summary(t19$d10.19[t19tr])
  t19[t19tr,c('Plot','Species','Num','d10.19')]
  
  # newr <- table(tAll$fsCat[r9],tAll$Type.19[r9])
  # print(newr)
  # s18 <- table(tAllm$fsCat[which(tAllm$Live.18==1)],tAllm$Type.19[which(tAllm$Live.18==1)])
  # #newr/s18
  # 
  # print('saplings by species and plot for unburned;low+medium;high severity')
  # r9ubsa <- intersect(r9,which(tAll$fsCat %in% c(0) & tAll$Type.19=='SA'))
  # print(table(tAll$Species[r9ubsa],tAll$Plot[r9ubsa]))
  # 
  # r9lmssa <- intersect(r9,which(tAll$fsCat %in% c(1,2) & tAll$Type.19=='SA'))
  # print(table(tAll$Species[r9lmssa],tAll$Plot[r9lmssa]))
  # 
  # r9hssa <- intersect(r9,which(tAll$fsCat==3 & tAll$Type.19=='SA'))
  # print(table(tAll$Species[r9hssa],tAll$Plot[r9hssa]))
  # 
  # What was the species composition of plots where Ceanothus germinated
  (cp <- unique(t19$Plot[which(t19$Species=='CEACUN')]))
  table(tAll$Species[which(tAll$Plot==cp[[1]])])
  table(tAll$Species[which(tAll$Plot==cp[[2]])])
  
  (cp <- unique(t19$Plot[which(t19$Species=='AMOCAL')]))
  table(tAll$Species[which(tAll$Plot %in% cp)])
}

write.csv(tAll,'data/tAll30.csv')

#### END HERE FOR NOW



# Pick a species, and run models in script 31

live.only <- F
yvar='gCrown.18'
#  run interactively in script 31
#P50r <- fitFatesMod(data=d,sp=spSel,fs='low-medium',logt=T)


# d <- tdat[[2]]
# table(d$Species)
# P50r <- fitARCMANmod()
# P50r
# (P50res <- rbind(P50res,P50r))

# sprouters
fst12c
dim(tAllm)

if (FALSE) {
  (spSel <- spAtt$Species[which(spAtt$Resprout=='Y' & spAtt$Shrub.Tree=='S')])
  spSel <- spSel[-which(spSel=='BACPIL')]
  spName <- 'RespShrubs'
  d <- barplotSprouterSpecies(spSel,FALSE,T,ss.name=spName)
  table(d$Species)
  fs=c('drop-high','low-medium') #'all','low-medium','drop-high'
  logt=T
  yvar='gCrown.18'
  # run interactively - fitFatesMod()
  
  d <- barplotSprouterSpecies(spSel,FALSE,T,ss.name=spName)
  yvar='DR.18'
  # run interactively - fitFatesMod()
  
  d <- barplotSprouterSpecies(spSel,FALSE,T,ss.name=spName)
  yvar='Live.18'
  # run interactively - fitFatesMod()
  
  d <- barplotSprouterSpecies(spSel,FALSE,T,ss.name=spName)
  # run multinomial interactively - fitMultiNomMod
  
  
  (spSel <- spAtt$Species[which(spAtt$Resprout=='Y' & spAtt$Shrub.Tree=='T')])
  spName <- 'RespTrees'
  d <- barplotSprouterSpecies(spSel,FALSE,T,ss.name=spName)
  table(d$Species)
  fs=c('low-medium') #'all','low-medium','drop-high'
  logt=T
  yvar='gCrown.18'
  # run interactively - fitFatesMod()
  
  d <- barplotSprouterSpecies(spSel,FALSE,T,ss.name=spName)
  # run multinomial interactively - fitMultiNomMod
}

# cycle through species here
spSel <- 'AMOCAL'
spSel <- 'UMBCAL'
spSel <- 'HETARB'
spSel <- 'QUEGAR'
spSel <- 'QUEAGR'

spName <- spSel
d <- barplotSprouterSpecies(spSel,FALSE,T)
table(d$Species)
fs='low-medium' #'all','low-medium','drop-high'
logt=T
live.only=F
spName <- spSel
# run interactively - fitFates2StepsMod()

yvar='Live.18'
live.only <- F
# run interactively - fitFatesMod()

d <- barplotSprouterSpecies(spSel,FALSE)
yvar='gCrown.18'
live.only <- T
# run interactively - fitFatesMod()

d <- barplotSprouterSpecies(spSel,FALSE)
yvar='DR.18'
# run interactively - fitFatesMod()

d <- barplotSprouterSpecies(spSel,FALSE)
# run multinomial interactively - fitMultiNomMod







#### UPDATE TO NEW FUNCTION AFTER HERE
spSel <- 'QUEGAR'
d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'gCrown'
# run model interactively

d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'DR'
# run model interactively

d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'Live'
# run model interactively


spSel <- 'HETARB'
d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
table(d$fsCat)
yvar <- 'gCrown'
# run model interactively

d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'DR'
# run model interactively

d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'Live'
# run model interactively

# Make a good set of bargraphs for selected species: PSEMEN, QUEAGR, UMBCAL, HETARB, QUEGAR

op=par(mfrow=c(5,4))
spSel <- 'PSEMEN'
tdat <- barplotOneNonSprouter(d,spSel,F,F,T)
spSel <- 'UMBCAL'
d <- barplotSprouterSpecies(spSel,F,F,T)
spSel <- 'QUEAGR'
d <- barplotSprouterSpecies(spSel,F,F,T)
spSel <- 'QUEGAR'
d <- barplotSprouterSpecies(spSel,F,F,T)
spSel <- 'HETARB'
d <- barplotSprouterSpecies(spSel,F,F,T)
par(op)



#### some continuous regressions
# Check histograms of sample sizes to decide which species can be analyzed - for d10 and sapling height
basalDiameterHistograms()

# just run logistic regressions on tree
basalDiameterLogisticRegressions()



##### TO HERE
# MOdels for ARCMAN, AMOCAL, ARBMEN - won't be using them
spSel <- 'AMOCAL'
d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar='gCrown'
# run fitAMOCALmod interactively
#(P50r <- fitAMOCALmod(yvar='gCrown')) # gCrown or DR for topkill

d <- barplotSprouterSpecies(spSel,FALSE)
yvar='DR'
# run fitAMOCALmod interactively
#(P50r <- fitAMOCALmod(yvar='DR')) # gCrown or DR for topkill
#(P50res <- rbind(P50res,P50r))

spSel <- 'ARBMEN'
d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'gCrown'

d <- barplotSprouterSpecies(spSel,FALSE)
head(d)
table(d$Species)
yvar <- 'DR'
#(P50r <- fitARBMENmod('gCrown')) # gCrown or DR for topkill
#(P50res <- rbind(P50res,P50r))

# barplots for common EHRO, WO, and resprouting shrubs
tdat <- barplotSprouters()
# functions don't converge for WHTO or RSHR, only for EHRO
d <- tdat[[2]]
fitSpeciesModelsContSize.gCrown()
fitSpeciesModelsContSize.Live()

# model fitting - resprouters
# need to go into the functions script and run this interactively, for model selection
# no results written to file for summarizing
fitModelsDiscreteFates()




# plot resilience??

#### END HERE 3/28/24

# # binomial  regression of 'live' and 'green crown' - won't try topkill as quadratic didn't work - in progress
# 
# tAll$gCrown.18 <- apply(tAll[,c('LN.18','LR.18')],1,max,na.rm=T)
# table(tAll$gCrown.18)
# 
# spsel <- 'PSEMEN'
# yvalname <- 'gCrown.18'
# 
# tAlls <- tAll[which(tAll$Species==spsel),]
# tAlls$yval <- tAlls[,yvalname]
# 
# fit5 <- glm(yval~ld10+ld10.2+fsCat+northness+SpCd14,data=tAlls,family='binomial')
# BIC(fit5)
# coefficients(fit5)
# 
# 
# 
# # Summary across types #GHYTIGYFYGIGH - all fates in a table for TR, SA and overall- all species included
# # I added spcode and species name to fst12, which broke code that used column numbers. So I've replaced them with column names, here and elsewhere below.
# tree.sum <- apply(fst12[fst12$Type=='TR',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
# tree.sum
# (tree.sum)/(tree.sum[1])
# 
# sap.sum <- apply(fst12[fst12$Type=='SA',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
# sap.sum
# (sap.sum)/(sap.sum[1])
# 
# all.sum <- tree.sum+sap.sum
# all.sum
# all.sum/all.sum[1]
# #ts.sum <- apply(fst12[fst12$Type=='TS',-c(1:2)],2,sum,na.rm=T)
# #(ts.sum-sap.sum[6])/(ts.sum[1]-ts.sum[6])
# 
# fate.sum <- apply(fst12[,-c(1:2)],2,sum)
# (fate.sum-fate.sum[6])/(fate.sum[1]-fate.sum[6])
# 
# SArows <- which(fst12$Type=='SA')
# TRrows <- which(fst12$Type=='TR')
# #TSrows <- which(fst12$Type=='TS') 
# 
# fst12$percSurv <- 1 - fst12$N18.DN/fst12$N17
# fst12$percSurv[fst12$N17==0] <- NA
# head(fst12)
# 
# # Add SpCd14 variable, with "Other" for everything that isn't the 9 primary
# tAll$SpCd14 <- tAll$Species.18
# 
# ## now run by species for species with lots of data
# spN <- table(tAll$SpCd14)
# spN[order(spN)]
# (spA <- names(spN)[which(spN>=50)])
# 
# # recode uncommon (<50) to 'Other
# tAll$SpCd14[which(!tAll$SpCd14 %in% spA & !is.na(tAll$Species.18))] <- 'OTHER'
# table(tAll$SpCd14)
# spA <- c(spA,'OTHER')
# 
# # Rerun outcomes table with 13 species + Other
# fst12a <- data.frame(SpCd14=rep(spA,each=2),Type=rep(c('SA','TR'),length(spA)),N17=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
# dim(fst12a)
# head(fst12a)                  
# tail(fst12a)
# 
# i=5
# for (i in 1:nrow(fst12a))
# {
#   sp <- fst12a$SpCd14[i]
#   ty <- fst12a$Type[i]
#   d <- tAll[which(tAll$SpCd14==sp & tAll$Type.17==ty),]
#   
#   fst12a$N17[i] <- sum(d$Live.17,na.rm=T)
#   
#   ## The next three lines are all equivalent - just using third one
#   #fst12a$N18.DN[i] <- length(which(d$fate.18=='DN'))
#   #fst12a$N18.DN[i] <- length(which(d$DN.18=='1'))
#   fst12a$N18.DN[i] <- sum(d$DN.18,na.rm = T)
#   
#   fst12a$N18.DR[i] <- sum(d$DR.18,na.rm = T)
#   fst12a$N18.LN[i] <- sum(d$LN.18,na.rm = T)
#   fst12a$N18.LR[i] <- sum(d$LR.18,na.rm = T)
#   miss <- which(d$Live.13==1 & is.na(d$DN.18)==1)
#   if (length(miss)>0) for (j in 1:length(miss)) print(d[miss[j],c('Plot.13','Num')])
#   fst12a$nMissing <- fst12a$N17-(fst12a$N18.DN+fst12a$N18.DR+fst12a$N18.LN+fst12a$N18.LR)
# }
# 
# fst12a
# head(fst12a)
# tail(fst12a)
# sum(fst12a$nMissing) #fixed 1330 dups we were previously ignoring- should now be zero
# 
# # Summary across types #GHYTIGYFYGIGH - all fates in a table for TR, SA and overall- all species included
# tree.sum <- apply(fst12a[fst12a$Type=='TR',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
# tree.sum
# (tree.sum)/(tree.sum[1])
# 
# sap.sum <- apply(fst12a[fst12a$Type=='SA',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
# sap.sum
# (sap.sum)/(sap.sum[1])
# 
# all.sum <- tree.sum+sap.sum
# all.sum
# all.sum/all.sum[1]
# 
# fate.sum <- apply(fst12a[,-c(1:2)],2,sum)
# (fate.sum-fate.sum[6])/(fate.sum[1]-fate.sum[6])
# 
# SArows <- which(fst12a$Type=='SA')
# TRrows <- which(fst12a$Type=='TR')
# #TSrows <- which(fst12a$Type=='TS') 
# 
# fst12a$percSurv <- 1 - fst12a$N18.DN/fst12a$N17
# fst12a$percSurv[fst12a$N17==0] <- NA
# head(fst12a)
# 
# # Table B - convert outcomes to percentages
# fst12a[fst12a$Type=='TR',]
# fst12a[fst12a$Type=='SA',]
# 
# ## choose fire severity metric
# fsmet <- 'Tubbs.MTBS.RDNBR.30'
# names(fs)
# summary(fs[,fsmet])
# hist(fs[,fsmet])
# sort(fs[,fsmet])
# 
# f2t <- match(tAll$Plot,fs$Plot)
# head(f2t)
# tail(f2t)
# tAll$FireSev <- fs[f2t,fsmet]
# dim(tAll)
# tail(tAll)
# 
# # RDNBR fire severity levels
# # Unburned
# # Low <170, but burned
# # Medium 170-700
# # High: >700
# # unburned: 1308, 1309, 1311, 1312, 1327, 1344, 1347
# 
# # Create a discrete FireSev variable
# summary(tAll$FireSev)
# hist(tAll$FireSev,breaks=c(-100,0,50,100,150,200,300,400,500,600,700,800,1000))
# fs.breaks <- c(-100,135,430,1000) # Based on Parks et al. 2014, but not using intermediate split at 304
# tAll$fsCat <- cut(tAll$FireSev,fs.breaks)
# table(tAll$fsCat)
# tAll$fsLevel <- as.numeric(tAll$fsCat)
# table(tAll$fsLevel)
# 
# # manually code unburned as 0
# unburned.plots <- c('PPW1308','PPW1309','PPW1311','PPW1312','PPW1327','PPW1344','PPW1347')
# tAll$fsLevel[which(tAll$Plot %in% unburned.plots)] <- 0
# table(tAll$fsLevel)
# tAll$fsCat <- as.factor(tAll$fsLevel)
# table(tAll$fsCat)
# str(tAll$fsCat)
# 
# # plots experiencing each fire severity level, and how many N13 individuals in each
# table(tAll$Plot.17[which(tAll$fsLevel==0)])
# table(tAll$Plot.17[which(tAll$fsLevel==1)])
# table(tAll$Plot.17[which(tAll$fsLevel==2)])
# table(tAll$Plot.17[which(tAll$fsLevel==3)])
# 
# # check taht only one fire severity per plot
# table(tAll$Plot,tAll$fsLevel)
# 
# # # Initial examination of dbh
# # dim(tAll)
# # length(which(tAll$Type.18=='TS')) #kuyihuiughiu
# # nodbh <- which(is.na(tAll$dbh.13))
# # length(nodbh)
# # head(tAll[nodbh,])
# # summary(tAll$dbh.13)
# 
# ## Use d10 for analysis - remember for SA is original basal data and for TR is calculated from dbh
# # D10 = DBH.cm * 1.176 + 1.070
# 
# # if want to change and use size from a different year, change here. Then from here on ld10 is generic
# hist(tAll$d10.17)
# tAll$ld10 <- log10(tAll$d10.17)
# tAll$ld10[which(!is.finite(tAll$ld10))] <- NA
# hist(tAll$ld10)
# summary(tAll$ld10,useNA='always')
# # smallest tree now had d10 = 1*1.176+1.07 = 2.246 for d10. log10 of this = 0.3514
# 
# # We've now changed to using 2018 size data unless missing, and then using 2013 - this comment and 5 lines of code below are from prior analysis. Left here as reminder about how this influenced the U-shaped multinomials. 
# 
# # now move diameters forward from 2013 (or 2018, if we use 2013 above), if it's missing in 2018. This can be commented out so that we only use data from one year or the other, and don't mix. I just tried this, for LN.18 analysis - it gets rid of the U-shaped result. So the result is coming from mixing data from plants censuses only in 2013 with those added in 2018. Now we have to figure out why!
# # which(is.na(tAll$ld10) & !is.na(tAll$d10.13))
# # tAll$ld10[which(is.na(tAll$ld10))] <- log10(tAll$d10.13[which(is.na(tAll$ld10))])
# # hist(tAll$ld10)
# 
# # create types to switch between types
# types <- c('SA','TR')
# op=par(mfrow=c(1,2))
# for (i in 1:2) {
#   ty <- types[i]
#   print(hist(tAll$ld10[tAll$Type.18==ty],main=paste('Type',ty)))
# }
# par(op)
# 
# # create squared variable for quadratic analysis
# tAll$ld10.2 <- tAll$ld10^2
# 
# # check on fates (for all 6945-703 indvs)
# table(tAll$Live.18, tAll$fate.18,useNA='always')
# 
# # gCrown not properly coded, fixed here using fates
# table(tAll$gCrown.18, tAll$fate.18,useNA='always')
# tAll$gCrown.18[which(tAll$fate.18 %in% c('LN','LR'))] <- 1
# tAll$gCrown.18[which(tAll$fate.18 %in% c('DN','DR'))] <- 0
# table(tAll$gCrown.18, tAll$fate.18,useNA='always')
# 
# table(tAll$bSprout.18,tAll$fate.18,useNA='always')
# 
# # 703 plants have fate=NA in 2018, let's have a look
# # these appear to be points added in 2019 and 2020 - that many? - Yes, especially lots of CEACUN!
# NA704 <- which(is.na(tAll$fate.18))
# tAll[NA704[3],]
# 
# # now try saplings and trees separately
# op <- par(mfrow=c(1,2))
# for (i in 1:2) plot(Live.18~ld10,data=tAll[tAll$Type.18==types[i],],main=paste('Type',types[i]))
# par(op)
# 
# # or plot together
# plot(Live.18~ld10,data=tAll,main=paste('All',types[i]))
# 
# # check fate values (for 6241 indvs- excluding the 703 NA's)
# table(tAll$Resprout.18,tAll$fate.18)
# 
# #### SURVIVAL ANALYSIS - first cut, size only!!
# types <- c('TR','SA')
# 
# # adjust values here to subset data, for TR and/or SA and fire severity level and species (if needed)
# table(tAll$Type.18,useNA='always')
# table(tAll$fsLevel,tAll$Plot,useNA='always')
# 
# write.csv(tAll,'data/tAll-analyzeData-update.csv')
# 
