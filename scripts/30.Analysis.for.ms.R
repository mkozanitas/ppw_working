## This is final analysis script with code in order required to reproduce results in Results section of manuscript

rm(list=ls())
#require(MuMln)
require(lme4)
require(glmmTMB)
require(Ternary)
require(vegan)
require(RCurl)
require(nnet)
require(brms)
#require(cmdstanr)

source('scripts/31.functionsForAnalysis.R')
op.reset <- par(mfcol=c(1,1),mar=c(5,5,3,3))
fsCols <- c('blue','brown','orange','red')

fs <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodFireSeverity/master/data/FSextract/vegplots-54-20m-FS.csv")
head(fs)
fs <- loadFireSeverity('Tubbs.MTBS.RDNBR.30',verbose=T)
head(fs[,c('Plot','fsvar','fsCat')])
table(fs$fsCat)

# Read in species codes - in this script, called 'Species'
spAtt <- read.csv('data/all-spp-attributes.csv',row.names = 1)
head(spAtt)
tail(spAtt)

# input merged dataframe
tAll <- read.csv('data/tAll.csv',as.is=T,row.names=1)
head(tAll)
dim(tAll)

# convert QSpecies# convert QUEBEGA to QUEBER, QUEDOGA to QUEDOU, QUEdec to QUEDOU, and remove unknowns
tAll <- convertHybrids()
table(tAll$Species)

# remove Type.13 = 'TS'; a few show up in 2018 - need to check
tAll <- tAll[-which(tAll$Type.17=='TS'),]
dim(tAll)

# add fire severity to tAll
tAll$fsCat <- fs$fsCat[match(tAll$Plot,fs$Plot)]
table(tAll$fsCat,useNA='always')

# archive tAll in tAll.arch and then reduce tAll to allow for easier examination during analysis
tAll.archive <- tAll

# use this to restore and recreate tAll
tAll <- tAll.archive
tAll <- tAll[,c('Num','iNum.13','Point.13','Plot','fsCat','Species','Year.13','Year.17','Type.17','Live.17','DBH_cm.17','d10.17','SA.Height_cm.17','Year.18','Type.18','Live.18','fate.18','DN.18','DR.18','LN.18','LR.18','gCrown.18','DBH_cm.18','d10.18','Basal.Resprout.Count.18','Basal.Resprout.Height_cm.18','Year.19','Live.19','fate.19','Type.19','DBH_cm.19','d10.19','Basal.Resprout.Count.19','Basal.Resprout.Height_cm.19')]

# Introduction - lists number of stems sampled. Include 2018 plots as this is just a general statement about size of the network
sum(tAll$Live.17,na.rm=T)

# Number of live species prefire
length(table(tAll$Species[which(tAll$Live.17==1)]))
spAtt

# First result - ternary plot ordination of plots to describe veg types, intersected with fire severity
pt <- drawTernaryPlots()
if (all(pt$plot==fs$Plot)) pt$fsCat <- fs$fsCat else print('error')
head(pt)
table(pt$vt)
table(pt$fsCat)
(vtfs <- table(pt[,c('vt','fsCat')]))
write.csv(vtfs,'results/veg-type-fire-sev-table.csv')

# Second result - overall numbers
table(tAll$Type.17)
(allAb <- sum(table(tAll$Type.17[which(tAll$Live.17==1)])))
(spAbund <- sort(table(tAll$Species[which(tAll$Live.17==1)])))

#### Summary stats
(use.species <- spAtt$Species)

(common.species <- spAtt$Species[which(spAtt$Common=='Yes')])
(sumAb <- sum(spAbund[which(names(spAbund) %in% common.species)]))
sumAb/allAb

# MELINA: QUESTION
if (FALSE) {
  # how many 'new' saplings 
  newSap <- which(is.na(tAll$Year.13) & tAll$Year.18==2018 & tAll$Type.18=='SA')
  table(tAll$Plot[newSap])
  length(newSap)
  table(tAll$Year.13,tAll$Year.17,useNA='always')
  
  tAllm$Num[newSap]
  table(tAllm$Plot[newSap])
  table(tAllm$fate.18[newSap],tAllm$fsCat[newSap])
  
  # MELINA: Remind me about these - they were either tagged but data failed to record in 2013, or newly recruited from 2013 to 2017 and inferred in 2017 data. They are detected here by NA in Year.17 which I think means 
}

# summary of species types
table(spAtt$Shrub.Tree)

# summary of fates
(ftab <- table(tAll$fate.18[which(tAll$Type.18!='TS')]))
(allAb18 <- sum(table(tAll$fate.18[which(tAll$Type.18!='TS')])))
ftab/allAb18
sum((ftab/allAb18)[3:4])

# create fst dataframe - FateSummaryTable for time 1 -> 2 (2017 and 2018)
fst12 <- calcFatesTableBySpecies()
head(fst12)
tail(fst12)

# reduce fst12 table to common species plus other shrubs and trees each tallied by SA and TR
fst12c <- reduce_fst12()
fst12c

write.csv(fst12c,'results/summary-table.csv')

#### make bar plots for different subgroups and individual species
# tAllm has all the living plants from 2017, with their 2018 fates
# doesn't have new plants added in 2019
tAllm <- prepareForBarPlots()
dim(tAll)
dim(tAllm)

# barplot for PSEMEN and ARCMAN, first to see here, and then to pdf
tdat <- barplotNonSprouters()

# model fitting - PSEMEN
  d <- tdat[[1]]
  head(d)
  table(d$Species)
  # Model comparison in here - can skip 
  {
    # obsolete models
   {  
  fit2 <- glm(gCrown.18~SizeCat + fsCat + SizeCat*fsCat,data=temp,family='binomial')
  summary(fit2)
  AIC(fit2)
  drop1(fit2)
  
  fit1 <- glm(gCrown.18~SizeCat + fsCat ,data=temp,family='binomial')
  summary(fit1)
  AIC(fit1)
  drop1(fit1)
  
  fit3 <- glmer(gCrown.18~d10.17 + fsCat + d10.17*fsCat + (1|Plot),data=temp,family='binomial')
  summary(fit3)
  AIC(fit3)
  
  fit3 <- lme4::glmer(gCrown.18~d10.17 + fsCat + d10.17*fsCat + (1|Plot),data=temp,family='binomial')
  summary(fit3)
  AIC(fit3)
  
  fit3t <- glmmTMB(gCrown.18~d10.17 + fsCat + d10.17*fsCat + (1|Plot),data=temp,family='binomial')
  summary(fit3t)
  AIC(fit3t)
  
  fit4 <- lme4::glmer(gCrown.18~d10.17 + fsCat + d10.17*fsCat + (1|Plot) + (1|iNum.13),data=temp,family='binomial')
  summary(fit4)
  
  fit4 <- glmmTMB(gCrown.18~d10.17 + fsCat + d10.17*fsCat + (1|Plot) + (1|iNum.13),data=temp,family='binomial')
  summary(fit4)
  AIC(fit4)
}
  fit5d <- glmmTMB(gCrown.18~factor(TSizeCat)*factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5d)
  AIC(fit5d)
  
  fit5 <- glmmTMB(gCrown.18~d10.17*factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5)
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(gCrown.18~d10.17 + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5_for_lrt)
  anova(fit5, fit5_for_lrt, "LRT")
  AIC(fit5_for_lrt)
  
  fit5_for_lrt2 <- glmmTMB(gCrown.18~factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5_for_lrt2)
  anova(fit5, fit5_for_lrt2, "LRT")
  AIC(fit5_for_lrt2)
  
  fit5_for_lrt3 <- glmmTMB(gCrown.18~d10.17 + factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5_for_lrt3)
  anova(fit5, fit5_for_lrt3, "LRT")
  AIC(fit5_for_lrt3)
  
  # fit6 <- brm(gCrown.18 ~
  #               d10.17*mo(as.integer(fsCat)) + # mo() treats the fire severity effect as monotonic
  #               (1|Plot)+
  #               (1|iNum.13),
  #             data=temp,
  #             family="Bernoulli",
  #             chains=2, # run two Stan sampler chains
  #             cores=2, # use two cores, so the two chains can run in parallel
  #             seed=238, # for reprodocibility
  #             #backend = "cmdstanr"  # this makes the output easier to read on some operating systems
  #             )
  # summary(fit6) # positive size effect, negative severity effect, negative interaction.
  # 
  # # quick viz of main effects and interaction
  # conditional_effects(fit6)
}
  # fit5 is best model - both terms supported in best model, including interaction, even if not individually significant in summary statement
  fit5 <- glmmTMB(gCrown.18~d10.17*factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
summary(fit5)

# model fitting - ARCMAN
{
  d <- tdat[[2]]
  #obsolete models
  {
    fit2 <- glm(gCrown.18~SizeCat + fsCat + SizeCat*fsCat,data=temp,family='binomial')
  AIC(fit2)
  drop1(fit2)
  fit1 <- glm(gCrown.18~SizeCat + fsCat,data=temp,family='binomial')
  AIC(fit1)
  drop1(fit1)
  fit0 <- glm(gCrown.18~fsCat ,data=temp,family='binomial')
  AIC(fit0)
  }
  fit5d <- glmmTMB(gCrown.18~factor(SSizeCat)*factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5d)
  AIC(fit5d)
  
  fit5 <- glmmTMB(gCrown.18~d10.17*factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5)
  AIC(fit5)
  
  fit5_for_lrt <- glmmTMB(gCrown.18~d10.17 + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5_for_lrt)
  anova(fit5, fit5_for_lrt, "LRT")
  AIC(fit5_for_lrt)
  
  fit5_for_lrt2 <- glmmTMB(gCrown.18~factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5_for_lrt2)
  anova(fit5, fit5_for_lrt2, "LRT")
  AIC(fit5_for_lrt2)
  
  fit5_for_lrt3 <- glmmTMB(gCrown.18~d10.17 + factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
  summary(fit5_for_lrt3)
  anova(fit5, fit5_for_lrt3, "LRT")
  AIC(fit5_for_lrt3)
}
# Best fit model or ARCMAN - no interaction of size and FS - summary can't calculate significance values, but model comparison worked to identify as best fit
fit5_for_lrt3 <- glmmTMB(gCrown.18~d10.17 + factor(fsCat) + (1|Plot) + (1|iNum.13), data=d, family='binomial')
summary(fit5_for_lrt3)

##### TO HERE
# barplots for common EHRO, WO, and resprouting shrubs
tdat <- barplotSprouters()

# model fitting - resprouters
# need to go into the functions script and run this interactively, for model selection
# no results written to file for summarizing
fitModelsDiscreteFates()

## what about delayed mortality
delayedMortality()

#### some continuous regressions
# Check histograms of sample sizes to decide which species can be analyzed - for d10 and sapling height
basalDiameterHistograms()

# just run logistic regressions on tree
basalDiameterLogisticRegressions()

## resprout height growth
basalResprouts()

# 2019 recruits
length(which(tAllm$Live.18==1))
length(which(tAllm$Live.18==1 & tAllm$Type.18=='SA'))
new2019recruits()


# plot resilience??

#### END HERE 3/28/26

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
#   temp <- tAll[which(tAll$SpCd14==sp & tAll$Type.17==ty),]
#   
#   fst12a$N17[i] <- sum(temp$Live.17,na.rm=T)
#   
#   ## The next three lines are all equivalent - just using third one
#   #fst12a$N18.DN[i] <- length(which(temp$fate.18=='DN'))
#   #fst12a$N18.DN[i] <- length(which(temp$DN.18=='1'))
#   fst12a$N18.DN[i] <- sum(temp$DN.18,na.rm = T)
#   
#   fst12a$N18.DR[i] <- sum(temp$DR.18,na.rm = T)
#   fst12a$N18.LN[i] <- sum(temp$LN.18,na.rm = T)
#   fst12a$N18.LR[i] <- sum(temp$LR.18,na.rm = T)
#   miss <- which(temp$Live.13==1 & is.na(temp$DN.18)==1)
#   if (length(miss)>0) for (j in 1:length(miss)) print(temp[miss[j],c('Plot.13','Num')])
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
