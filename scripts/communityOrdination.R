#analyzeData-tAll is to compare 2013, 2018, 2019, and 2020 for delayed mortality and other changes in state

# RUN prepareData and examineData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R

rm(list=ls())
library(lme4)
library(glmmTMB)
require(Ternary)

source('scripts/31.functionsForAnalysis.R')
fsCols <- c('blue','brown','orange','red')

# Functions to convert basal area to diameter and log-diameter, and back
lba2d <- function(x) 2*sqrt((10^x)/pi)
ba2d <- function(x) 2*sqrt((x)/pi)
d2ba <- function(x) pi*(x/2)^2
d2lba <- function(x) log10(pi*(x/2)^2)

# Load fire severity data
library("RCurl")
fs <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodFireSeverity/master/data/FSextract/vegplots-54-20m-FS.csv")
head(fs)
tail(fs)

fsmet <- 'Tubbs.MTBS.RDNBR.30'
names(fs)
fs$Plot
summary(fs[,fsmet])
hist(fs[,fsmet])
sort(fs[,fsmet])

# RDNBR fire severity levels
# Unburned
# Low <170, but burned
# Medium 170-700
# High: >700
# unburned: 1308, 1309, 1311, 1312, 1327, 1344, 1347

# Create a discrete FireSev variable
hist(fs[,fsmet],breaks=c(-100,0,50,100,150,200,300,400,500,600,700,800,1000))
fs.breaks <- c(-100,135,430,1000) # Based on Parks et al. 2014, but not using intermediate split at 304
fs$fsCat <- cut(fs[,fsmet],fs.breaks)
table(fs$fsCat)
fs$fsLevel <- as.numeric(fs$fsCat)
table(fs$fsLevel)

# manually code unburned as 0
unburned.plots <- c('PPW1308','PPW1309','PPW1311','PPW1312','PPW1327','PPW1344','PPW1347')
fs$fsLevel[which(fs$Plot %in% unburned.plots)] <- 0
table(fs$fsLevel)
fs$fsCat <- as.factor(fs$fsLevel)
table(fs$fsCat)
str(fs$fsCat)
####

# Read in species codes - in this script, called 'Species'
spNames <- read.csv('data/all-spp-names.csv',row.names = 1)
head(spNames)
names(spNames)[1] <- 'Species'

# input merged dataframe
tAll <- read.csv('data/tAll.csv',as.is=T,row.names=1)
head(tAll)
dim(tAll)
dim(tAll)
names(tAll)

tAlls <- tAll[,c('Num','Plot','Species','Type.17','d10.17','fate.18','Type.18','d10.18','Basal.Resprout.Count.18','fate.19','Type.19','d10.19','DBH_cm.19','Basal.Resprout.Count.19','fate.20','Type.20','d10.20','Basal.Resprout.Count.20')]
(brc <- grep('Basal.Resprout.Count',names(tAlls)))
names(tAlls)[brc] <- c('BRC.18','BRC.19','BRC.20')
# 
(tAlls[which(tAlls$Plot=='PPW1302' & tAlls$fate.18=='DR'),])

table(tAll$fate.18,tAll$gCrown.18)

# ternary plot of the three major groups
spAtt <- read.csv('data/all-spp-attributes.csv')
head(spAtt)
s2sg <- match(tAlls$Species,spAtt$OrigSpecies)
isna <- which(is.na(s2sg))
tAlls$Species[isna]
tAlls$SGroup <- spAtt$SGroup[s2sg]
tAlls$BA.17 <- d2ba(tAlls$d10.17)

pt <- drawTernaryPlots()

SGBA <- tapply(tAlls$BA.17,list(tAlls$Plot,tAlls$SGroup),sum,na.rm=T)
SGBA[is.na(SGBA)] <- 0
sumBA <- apply(SGBA,1,sum)
pSGBA <- SGBA/sumBA
pSGBA3 <- pSGBA
pSGBA3[,4] <- pSGBA3[,3] + pSGBA3[,4]
(pSGBA3 <- round(100*pSGBA3[,-3],1))

# classify plots
pt <- data.frame(plot=row.names(pSGBA),vt=NA)
pt$vt[which(pSGBA3[,1]>=50)] <- 'CON'
pt$vt[which(pSGBA3[,2]>=50)] <- 'MH'
pt$vt[which(pSGBA3[,3]>=50)] <- 'WO'
pt$vt[which(is.na(pt$vt))] <- 'Mix3'
pt
pt$fs <- fs$fsCat[match(pt$plot,fs$Plot)]
head(pt)
table(pt$vt,pt$fs)

library(vegan)
sgord <- decorana(pSGBA3)
plot(sgord)

## ## sum## ## sumpSGBA
# Structural ordination
# make a matrix where each plot has all the trees counted in structural categories, in each year. The goal is to see how the plot structure changes post-fire, using ordination
# categories: 4 sizes by 2 status = 8:
# saplings, small trees (1-10), medium trees (10-20), large trees (20+)
# intersected with post-fire status: green crown or topkill/resprout (dead trees are dropped); topkill/resprout are put into size categories based on pre-fire
# the final matrix should have 54 x 4 rows: plots by time periods, and 8 columns; 

# ps = plot structure
plots <- sort(unique(tAll$Plot.17))
years <- c('17','18','19','20')
ps <- data.frame(plot=rep(plots,each=4),year=rep(years,54),SAG=NA,STG=NA,MTG=NA,LTG=NA,SAR=NA,STR=NA,MTR=NA,LTR=NA,SAD=NA,STD=NA,MTD=NA,LTD=NA)
head(ps)


j=1
r <- 0
for (j in 1:length(plots))
{
  i=4
  for (i in 1:length(years)) 
  {
    r <- r+1
    #r <- (i-1)*length(plots) + j
    p <- plots[j]
    y <- years[i]
    tp <- paste('Plot.',y,sep='')
    dim(tAll)
    tx <- tAll[!is.na(tAll[,c(paste('Plot.',y,sep=''))]),]
    dim(tx)
    txx <- tx[tx[,tp]==p,]
    dim(txx)
    if (nrow(txx)>0) {
      if (i==1)
      {
        ps[r,'SAG'] <- length(which(txx$Type.17=='SA'))
        ps[r,'STG'] <- length(which(txx$Type.17=='TR' & txx$DBH_cm.17<=10))
        ps[r,'MTG'] <- length(which(txx$Type.17=='TR' & txx$DBH_cm.17>10  & txx$DBH_cm.17<=20))
        ps[r,'LTG'] <- length(which(txx$Type.17=='TR' & txx$DBH_cm.17>20))
        ps[r,'SAR'] <- 0
        ps[r,'STR'] <- 0
        ps[r,'MTR'] <- 0
        ps[r,'LTR'] <- 0
        ps[r,'SAD'] <- 0
        ps[r,'STD'] <- 0
        ps[r,'MTD'] <- 0
        ps[r,'LTD'] <- 0
      } else 
      {
        ty <- paste('Type.',y,sep='')
        td <- paste('DBH_cm.',y,sep='')
        tv <- paste('Topkill.',y,sep='')
        lv <- paste('Live.',y,sep='')
        t0 <- which(txx[,tv]==0 & txx[,lv]==1)
        t1 <- which(txx[,tv]==1 & txx[,lv]==1)
        t2 <- which(txx[,lv]==0)
        ps[r,'SAG'] <- length(intersect(
          t0,which(txx[,ty]=='SA')))
        ps[r,'STG'] <- length(intersect(
          t0,which(txx[,ty]=='TR' & txx[,td]<=10)))
        ps[r,'MTG'] <- length(intersect(
          t0,which(txx[,ty]=='TR' & txx[,td]>10  & txx[,td]<=20)))
        ps[r,'LTG'] <- length(intersect(
          t0,which(txx[,ty]=='TR' & txx[,td]>20)))
        ps[r,'SAR'] <- length(intersect(
          t1,which(txx[,ty]=='SA')))
        ps[r,'STR'] <- length(intersect(
          t1,which(txx[,ty]=='TR' & txx[,td]<=10)))
        ps[r,'MTR'] <- length(intersect(
          t1,which(txx[,ty]=='TR' & txx[,td]>10  & txx[,td]<=20)))
        ps[r,'LTR'] <- length(intersect(
          t1,which(txx[,ty]=='TR' & txx[,td]>20)))
        ps[r,'SAD'] <- length(intersect(
          t2,which(txx[,ty]=='SA')))
        ps[r,'STD'] <- length(intersect(
          t2,which(txx[,ty]=='TR' & txx[,td]<=10)))
        ps[r,'MTD'] <- length(intersect(
          t2,which(txx[,ty]=='TR' & txx[,td]>10  & txx[,td]<=20)))
        ps[r,'LTD'] <- length(intersect(
          t2,which(txx[,ty]=='TR' & txx[,td]>20)))
      }
    }
  }
}

ps$PY <- paste(ps$plot,ps$year,sep='.')
head(ps)
tail(ps)
ps[ps$plot=='PPW1346',]
write.csv(ps,'data/plot-structure-x-year.csv')

calcProp <- function(x) x/sum(x)
useProp <- T

library(vegan)
ps2 <- ps[ps$year %in% c(17:18),]
head(ps2[,3:14])
if (useProp) {
  for (i in 1:nrow(ps2)) ps2[i,3:14] <- calcProp(ps2[i,3:14])
  ps2[,3:14] <- round(ps2[,3:14]*100,0)
}
head(ps2[,3:14])

rownames(ps2) <- 1:nrow(ps2)
dim(ps2)
inclDead <- TRUE
if (inclDead) ords <- decorana(ps2[,3:14]) else ords <- decorana(ps2[,3:10])
plot(ords)
attributes(summary(ords))
ss <- summary(ords)$site.scores[,1:2]
head(ss)
tail(ss)
points(ss[ps2$year==17,],pch=19)
points(ss[ps2$year==18,],pch=19,col='red')
#points(ss[109:162,],pch=19,col='green')
dim(ss)
ss.init <- ss[ps2$year==17,]
head(ss.init)

i=1
psd <- data.frame(plot=plots,pn=1:54,vlen=NA,vang=NA)
for (i in 1:54)
{
  ii <- i*2-1
  r <- c(ii,ii+1)
  j=1
  for (j in 1:1) {
    x1 <- ss[c(r[j]),1]
    x2 <- ss[c(r[j+1]),1]
    y1 <- ss[c(r[j]),2]
    y2 <- ss[c(r[j+1]),2]
    psd$vlen[i] <- sqrt((x1-x2)^2 + (y1-y2)^2)
    if (psd$vlen[i]==0) psd$vang[i] <- 0 else 
      if ((x2-x1)==0) {
        if ((y2-y1)>0) psd$vang[i] <- 1 else psd$vang[i] <- -1
      }
    else 
      psd$vang[i] <- atan((y2-y1)/(x2-x1))
    lines(c(x1,x2),c(y1,y2))
  }
}
(psd)
plot(psd$vlen,psd$vang,type='n')
text(psd$vlen,psd$vang,labels=psd$pn)


psd$fsLevel <- fs$fsLevel
boxplot(psd$vlen~as.factor(psd$fsLevel))
boxplot(psd$vang~as.factor(psd$fsLevel))

# clusters based on structural ordination on absolute #s, not proportions
#cl <- read.csv('data/ordination/clusters.csv')
#head(cl)
#table(psd$fsLevel,cl$Cluster)

####
# How about species composition instead of structure
# do this two ways - separating saplings and trees, or combining into one count
splist <- sort(unique(tAll$Species))
(splist2 <- splist)
splist2[which(splist2=='QUEBEGA')] <- 'QUEBGA'
splist2[which(splist2=='QUEDOGA')] <- 'QUEDGA'
splist2[which(splist2=='SAL_SP')] <- 'SAL.SP'
splist2[which(splist2=='UNK')] <- 'UNKNOW'
splist2
tAll$Species2 <- splist2[match(tAll$Species,splist)]
table(tAll$Species2)

ty <- c('SA','TR')
st <- paste(rep(splist2,2),rep(ty,each=length(splist)),sep='_')
sp <- splist2
head(st)
length(st)

plots <- sort(unique(tAll$Plot.17))
years <- c('17','18','19','20')
pct <- data.frame(plot=rep(plots,each=4),year=rep(years,54))
for (i in 1:length(st)) pct <- cbind(pct,NA)
dim(pct)
names(pct)
names(pct)[-c(1:2)] <- st
head(pct)

pc <- data.frame(plot=rep(plots,each=4),year=rep(years,54))
for (i in 1:length(sp)) pc <- cbind(pc,NA)
dim(pc)
names(pc)
names(pc)[-c(1:2)] <- sp
head(pc)

r <- 0
j=1
for (j in 1:length(plots))
{
  i=1
  for (i in 1:length(years)) 
  {
    r <- r+1
    y <- years[i]
    dim(tAll)
    tx <- tAll[!is.na(tAll[,c(paste('Plot.',y,sep=''))]),]
    dim(tx)
    # r <- (i-1)*length(plots) + j
    p <- plots[j]
    tp <- paste('Plot.',y,sep='')
    tl <- paste('Live.',y,sep='')
    ty <- paste('Type.',y,sep='')
    txx <- tx[tx[,tp]==p & tx[,tl]==1,]
    dim(txx)
    if (nrow(txx)>0) {
      txx$st <- paste(txx$Species2,txx[,ty],sep='_')
      tst <- table(txx$st)
      x2x <- match(st,names(tst))
      t2t <- tst[x2x]
      t2t[which(is.na(t2t))] <- 0
      pct[r,-c(1:2)] <- t2t
      
      tst <- table(txx$Species2)
      x2x <- match(sp,names(tst))
      t2t <- tst[x2x]
      t2t[which(is.na(t2t))] <- 0
      pc[r,-c(1:2)] <- t2t
    }
  }
}
pct$PY <- paste(ps$plot,ps$year,sep='.')
head(pct)
tail(pct)
apply(pct[,-c(1,2,67)],2,sum,na.rm=T)

pc$PY <- paste(ps$plot,ps$year,sep='.')
head(pc)
tail(pc)
apply(pc[,-c(1,2,ncol(pc))],2,sum,na.rm=T)

# number of stems censused across year 1 and all four years!
sum(apply(pct[pct$year==17,-c(1,2,67)],2,sum,na.rm=T))
sum(apply(pct[,-c(1,2,67)],2,sum,na.rm=T))

dim(pct)
names(pct)

## IF DESIRED remove zero cols
# stem.sums <- apply(pct[pct$year<19,-c(1,2,67)],2,sum,na.rm=T)
# stem.sums
# 
# pct2 <- pct[pct$year<19,-(which(stem.sums==0)+2)]
# dim(pct2)
# apply(pct2[,-c(1,2,ncol(pct2))],2,sum,na.rm=T)

pct2 <- pct[pct$year<19,]
head(pct2)

pc2 <- pc[pc$year<19,]
head(pc2)

if (useProp) {
  for (i in 1:nrow(pct2)) pct2[i,-c(1,2,ncol(pct2))] <- calcProp(pct2[i,-c(1,2,ncol(pct2))])
  pct2[,-c(1,2,ncol(pct2))] <- round(pct2[,-c(1,2,ncol(pct2))]*100,0)
  head(pct2)
  
  for (i in 1:nrow(pc2)) pc2[i,-c(1,2,ncol(pc2))] <- calcProp(pc2[i,-c(1,2,ncol(pc2))])
  pc2[,-c(1,2,ncol(pc2))] <- round(pc2[,-c(1,2,ncol(pc2))]*100,0)
  head(pc2)
}

dim(pct2)
head(pct2)
ordct <- decorana(pct2[,-c(1,2,ncol(pct2))])
plot(ordct)
attributes(summary(ordct))
sptsc <- summary(ordct)$spec.scores[,1:2]
head(sptsc)
sptsc[order(sptsc[,1]),]

stc <- summary(ordct)$site.scores[,1:2]
points(stc[pct2$year==17,],pch=19)
points(stc[pct2$year==18,],pch=19,col='red')
#points(sc[pct2$year==19,],pch=19,col='green')
dim(stc)
stc.init <- stc[pct2$year==17,]

#calculate length of vectors connecting points
pctd <- data.frame(plot=plots,cvlen=NA,cvang=NA)

i=11
for (i in 1:54)
{
  ii <- i*2-1
  r <- c(ii,ii+1)
  j=1
  for (j in 1:1) {
    x1 <- stc[r[j],1]
    x2 <- stc[r[j+1],1]
    y1 <- stc[r[j],2]
    y2 <- stc[r[j+1],2]
    dd <- sqrt((x1-x2)^2+(y1-y2)^2)
    pctd$cvlen[i] <- dd
    if (pctd$cvlen[i]==0) pctd$cvang[i] <- 0 else 
      if ((x2-x1)==0) {
        if ((y2-y1)>0) pctd$cvang[i] <- 1 else pctd$cvang[i] <- -1
      }
    else 
      pctd$cvang[i] <- atan((y2-y1)/(x2-x1))
    lines(c(x1,x2),c(y1,y2))
  }
}
head(pctd)
tail(pctd)
pctd[order(pctd$cvlen),]

boxplot(pctd$cvlen~as.factor(psd$fsLevel))
plot(psd$vlen,pctd$cvlen)

# ordinate composition only
ordc <- decorana(pc2[,-c(1,2,ncol(pc2))])
plot(ordc)
attributes(summary(ordc))
spsc <- summary(ordc)$spec.scores[,1:2]
head(spsc)
spsc[order(spsc[,1]),]
plot(spsc,type='n')
text(spsc,labels = row.names(spsc))

plot(ordc)
sc <- summary(ordc)$site.scores[,1:2]
points(sc[pc2$year==17,],pch=19)
points(sc[pc2$year==18,],pch=19,col='red')
#points(sc[pct2$year==19,],pch=19,col='green')
dim(sc)
sc.init <- sc[pc2$year==17,]

#calculate length of vectors connecting points
pcd <- data.frame(plot=plots,cvlen=NA,cvang=NA)

i=11
for (i in 1:54)
{
  ii <- i*2-1
  r <- c(ii,ii+1)
  j=1
  for (j in 1:1) {
    x1 <- sc[r[j],1]
    x2 <- sc[r[j+1],1]
    y1 <- sc[r[j],2]
    y2 <- sc[r[j+1],2]
    dd <- sqrt((x1-x2)^2+(y1-y2)^2)
    pcd$cvlen[i] <- dd
    if (pcd$cvlen[i]==0) pcd$cvang[i] <- 0 else 
      if ((x2-x1)==0) {
      }
    else 
      pcd$cvang[i] <- atan((y2-y1)/(x2-x1))
    lines(c(x1,x2),c(y1,y2))
  }
}
head(pcd)
tail(pcd)
pcd[order(pcd$cvlen),]

boxplot(pcd$cvlen~as.factor(psd$fsLevel))
plot(psd$vlen,pcd$cvlen)

pairs(data.frame(pcd$cvlen,pctd$cvlen,psd$vlen))
head(psd)
psd$plot[which(psd$vlen>1.75 & pctd$cvlen<0.4)]

lvar <- pctd$cvlen
fit <- glm(lvar ~ factor(psd$fsLevel) + pt$vt)
summary(fit)

fit2 <- glm(lvar ~ factor(psd$fsLevel))
anova(fit,fit2,"LRT")

fit3 <- glm(lvar ~ pt$vt)
summary(fit3)
anova(fit,fit3,"LRT")

# order plots by compositional characteristics
head(pct2)
spat <- read.csv('data/all-spp-attributes.csv')
head(spat)
(spty <- names(pct2))
s2s <- match(substr(spty,1,6),spat$Species)
ST <- spat$Shrub.Tree[s2s]
DE <- spat$Dec.Ever[s2s]
RS <- spat$Resprout[s2s]

pcatt <- data.frame(plot=pct2$plot,year=pct2$year, Stot=NA,Ttot=NA,Dtot=NA,Etot=NA,Rtot=NA,Ntot=NA)
pcatt$Stot <- apply(pct2[,which(ST=='S')],1,sum)
pcatt$Ttot <- apply(pct2[,which(ST=='T')],1,sum)
pcatt$Dtot <- apply(pct2[,which(DE=='D')],1,sum)
pcatt$Etot <- apply(pct2[,which(DE=='E')],1,sum)
pcatt$Rtot <- apply(pct2[,which(RS=='Y')],1,sum)
pcatt$Ntot <- apply(pct2[,which(RS=='N')],1,sum)
head(pcatt)

pcatt$pT <- pcatt$Ttot/(pcatt$Ttot+pcatt$Stot)
pcatt$pE <- pcatt$Etot/(pcatt$Etot+pcatt$Dtot)
pcatt$pR <- pcatt$Rtot/(pcatt$Rtot+pcatt$Ntot)
head(pcatt)
tail(pcatt)

# ordination by functional attributes
pca2 <- pcatt[pcatt$year==17,]
row.names(pca2) <- 1:nrow(pca2)
pairs(pca2[,c('pT','pE','pR')])
orda <- decorana(pca2[,c('pT','pE','pR')])
plot(orda)

sa <- summary(orda)$site.scores[,1:2]
head(sa)
tail(sa)

## COMPILE ALL RESULTS
ar <- cbind(pca2[,c('pT','pE','pR')],sa,ss.init,stc.init,sc.init,pcd,pctd,psd)
dim(ar)
names(ar)
names(ar) <- c('pT','pE','pR','AttOrd1','AttOrd2','StrOrd1','StrOrd2','CStOrd1','CStOrd2','ComOrd1','ComOrd2','Plot','cvlen','cvang','plot','ctvlen','ctvang','plot2','pn','svlen','svang','FireSev')

plot(cvlen~pR,data=ar)

## Too many predictors!!!!
fit1 <- lm(cvlen~FireSev+AttOrd1+AttOrd2+StrOrd1+StrOrd2+CStOrd1+CStOrd2+ComOrd1+ComOrd2,data=ar)
summary(fit1)

# AttOrd2 positive effect - what are species loadings
sptsc[order(sptsc[,1]),]
spsc[order(spsc[,1]),]

fit2 <- lm(ctvlen~FireSev+AttOrd1+AttOrd2+StrOrd1+StrOrd2+CStOrd1+CStOrd2+ComOrd1+ComOrd2,data=ar)
summary(fit2)

fit3 <- lm(svlen~FireSev+AttOrd1+AttOrd2+StrOrd1+StrOrd2+CStOrd1+CStOrd2+ComOrd1+ComOrd2,data=ar)
summary(fit3)
plot(svlen~StrOrd1,data=ar)
summary(ords)$spec.scores

plot(AttOrd1~StrOrd1,data=ar)
plot(cvlen~AttOrd1,data=ar)
plot(ordct)


plot(svlen~StrOrd1,data=ar)
plot(ords)

plot(orda)

plot(svlen~fsLevel,data=ar)
plot(cvlen~fsLevel,data=ar)

########fit########################
# How many of the 2013s are TS - 250 - these are excluded from summary stats
length(which(tAll$Type.13=='TS'))
length(which(tAll$Type.17=='TS'))
length(which(tAll$Type.18=='TS'))
length(which(tAll$Type.19=='TS'))
length(which(tAll$Type.20=='TS'))

missStems <- which(is.na(tAll$Plot.18) & tAll$Type.19=='SA')
length(missStems)
head(tAll[missStems[1],])

table(tAll$Type.18,tAll$Type.19,useNA = 'always')

# initial growth analysis to identify potentially problematic size data
# suggest taking median eliminating most negative and very large outliers. Suggests median diameter growth of 0.2 cm in five years
tAll$ddbh.1318 <- tAll$dbh.18-tAll$dbh.13

plot(tAll$dbh.13,tAll$ddbh.1318)
abline(h=0)

tAll$absddbh.1318 <- abs(tAll$ddbh.1318)
hist(tAll$absddbh.1318)
median(tAll$ddbh.1318,na.rm=T)
tail(sort(tAll$absddbh.1318))

length(which(tAll$absddbh.1318>5))
write.csv(tAll[which(tAll$absddbh.1318>5),c('Num','Plot.18','ddbh.1318')],'data/growthquestions.csv')

## THIS CODE SECTION ANALYZES 2018 FATES ###

# ANALYZE BY SPECIES AND TYPE for 2018 post-fire fates
(use.species <- spNames$Species)

# hard code conversion of QUEBEGA to QUEBER
tAll$Species.13[which(tAll$Species.13=='QUEBEGA')] <- 'QUEBER'
tAll$Species.17[which(tAll$Species.17=='QUEBEGA')] <- 'QUEBER'
tAll$Species.18[which(tAll$Species.18=='QUEBEGA')] <- 'QUEBER'
tAll$Species.19[which(tAll$Species.19=='QUEBEGA')] <- 'QUEBER'
tAll$Species.20[which(tAll$Species.20=='QUEBEGA')] <- 'QUEBER'
tAll$Species[which(tAll$Species=='QUEBEGA')] <- 'QUEBER'

# create fst dataframe - FateSummaryTable for time 1 -> 2 (2017 and 2018)
fst12 <- data.frame(SpCode=rep(spNames$Species,each=2),Type=rep(c('SA','TR'),length(use.species)),N17=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
head(fst12)                  
tail(fst12)

i=5
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

fst12
head(fst12)
tail(fst12)
sum(fst12$nMissing) #fixed 1330 dups we were previously ignoring- should now be zero

# Summary across types #GHYTIGYFYGIGH - all fates in a table for TR, SA and overall- all species included
# I added spcode and species name to fst12, which broke code that used column numbers. So I've replaced them with column names, here and elsewhere below.
tree.sum <- apply(fst12[fst12$Type=='TR',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
tree.sum
(tree.sum)/(tree.sum[1])

sap.sum <- apply(fst12[fst12$Type=='SA',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
sap.sum
(sap.sum)/(sap.sum[1])

all.sum <- tree.sum+sap.sum
all.sum
all.sum/all.sum[1]
#ts.sum <- apply(fst12[fst12$Type=='TS',-c(1:2)],2,sum,na.rm=T)
#(ts.sum-sap.sum[6])/(ts.sum[1]-ts.sum[6])

fate.sum <- apply(fst12[,-c(1:2)],2,sum)
(fate.sum-fate.sum[6])/(fate.sum[1]-fate.sum[6])

SArows <- which(fst12$Type=='SA')
TRrows <- which(fst12$Type=='TR')
#TSrows <- which(fst12$Type=='TS') 

fst12$percSurv <- 1 - fst12$N18.DN/fst12$N17
fst12$percSurv[fst12$N17==0] <- NA
head(fst12)

# Add SpCd14 variable, with "Other" for everything that isn't the 9 primary
tAll$SpCd14 <- tAll$Species.18

## now run by species for species with lots of data
spN <- table(tAll$SpCd14)
spN[order(spN)]
(spA <- names(spN)[which(spN>=50)])

# recode uncommon (<50) to 'Other
tAll$SpCd14[which(!tAll$SpCd14 %in% spA & !is.na(tAll$Species.18))] <- 'OTHER'
table(tAll$SpCd14)
spA <- c(spA,'OTHER')

# Rerun outcomes table with 13 species + Other
fst12a <- data.frame(SpCd14=rep(spA,each=2),Type=rep(c('SA','TR'),length(spA)),N17=NA,N18.DN=NA,N18.DR=NA,N18.LN=NA,N18.LR=NA,nMissing=NA)
dim(fst12a)
head(fst12a)                  
tail(fst12a)

i=5
for (i in 1:nrow(fst12a))
{
  sp <- fst12a$SpCd14[i]
  ty <- fst12a$Type[i]
  temp <- tAll[which(tAll$SpCd14==sp & tAll$Type.17==ty),]
  
  fst12a$N17[i] <- sum(temp$Live.17,na.rm=T)
  
  ## The next three lines are all equivalent - just using third one
  #fst12a$N18.DN[i] <- length(which(temp$fate.18=='DN'))
  #fst12a$N18.DN[i] <- length(which(temp$DN.18=='1'))
  fst12a$N18.DN[i] <- sum(temp$DN.18,na.rm = T)
  
  fst12a$N18.DR[i] <- sum(temp$DR.18,na.rm = T)
  fst12a$N18.LN[i] <- sum(temp$LN.18,na.rm = T)
  fst12a$N18.LR[i] <- sum(temp$LR.18,na.rm = T)
  miss <- which(temp$Live.13==1 & is.na(temp$DN.18)==1)
  if (length(miss)>0) for (j in 1:length(miss)) print(temp[miss[j],c('Plot.13','Num')])
  fst12a$nMissing <- fst12a$N17-(fst12a$N18.DN+fst12a$N18.DR+fst12a$N18.LN+fst12a$N18.LR)
}

fst12a
head(fst12a)
tail(fst12a)
sum(fst12a$nMissing) #fixed 1330 dups we were previously ignoring- should now be zero

# Summary across types #GHYTIGYFYGIGH - all fates in a table for TR, SA and overall- all species included
tree.sum <- apply(fst12a[fst12a$Type=='TR',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
tree.sum
(tree.sum)/(tree.sum[1])

sap.sum <- apply(fst12a[fst12a$Type=='SA',c('N17','N18.DN','N18.DR','N18.LN','N18.LR')],2,sum)
sap.sum
(sap.sum)/(sap.sum[1])

all.sum <- tree.sum+sap.sum
all.sum
all.sum/all.sum[1]

fate.sum <- apply(fst12a[,-c(1:2)],2,sum)
(fate.sum-fate.sum[6])/(fate.sum[1]-fate.sum[6])

SArows <- which(fst12a$Type=='SA')
TRrows <- which(fst12a$Type=='TR')
#TSrows <- which(fst12a$Type=='TS') 

fst12a$percSurv <- 1 - fst12a$N18.DN/fst12a$N17
fst12a$percSurv[fst12a$N17==0] <- NA
head(fst12a)

# Table B - convert outcomes to percentages
fst12a[fst12a$Type=='TR',]
fst12a[fst12a$Type=='SA',]

## choose fire severity metric
fsmet <- 'Tubbs.MTBS.RDNBR.30'
names(fs)
summary(fs[,fsmet])
hist(fs[,fsmet])
sort(fs[,fsmet])

f2t <- match(tAll$Plot,fs$Plot)
head(f2t)
tail(f2t)
tAll$FireSev <- fs[f2t,fsmet]
dim(tAll)
tail(tAll)

# RDNBR fire severity levels
# Unburned
# Low <170, but burned
# Medium 170-700
# High: >700
# unburned: 1308, 1309, 1311, 1312, 1327, 1344, 1347

# Create a discrete FireSev variable
summary(tAll$FireSev)
hist(tAll$FireSev,breaks=c(-100,0,50,100,150,200,300,400,500,600,700,800,1000))
fs.breaks <- c(-100,135,430,1000) # Based on Parks et al. 2014, but not using intermediate split at 304
tAll$fsCat <- cut(tAll$FireSev,fs.breaks)
table(tAll$fsCat)
tAll$fsLevel <- as.numeric(tAll$fsCat)
table(tAll$fsLevel)

# manually code unburned as 0
unburned.plots <- c('PPW1308','PPW1309','PPW1311','PPW1312','PPW1327','PPW1344','PPW1347')
tAll$fsLevel[which(tAll$Plot %in% unburned.plots)] <- 0
table(tAll$fsLevel)
tAll$fsCat <- as.factor(tAll$fsLevel)
table(tAll$fsCat)
str(tAll$fsCat)

# plots experiencing each fire severity level, and how many N13 individuals in each
table(tAll$Plot.17[which(tAll$fsLevel==0)])
table(tAll$Plot.17[which(tAll$fsLevel==1)])
table(tAll$Plot.17[which(tAll$fsLevel==2)])
table(tAll$Plot.17[which(tAll$fsLevel==3)])

# check taht only one fire severity per plot
table(tAll$Plot,tAll$fsLevel)

# # Initial examination of dbh
# dim(tAll)
# length(which(tAll$Type.18=='TS')) #kuyihuiughiu
# nodbh <- which(is.na(tAll$dbh.13))
# length(nodbh)
# head(tAll[nodbh,])
# summary(tAll$dbh.13)

## Use d10 for analysis - remember for SA is original basal data and for TR is calculated from dbh
# D10 = DBH.cm * 1.176 + 1.070

# if want to change and use size from a different year, change here. Then from here on ld10 is generic
hist(tAll$d10.17)
tAll$ld10 <- log10(tAll$d10.17)
tAll$ld10[which(!is.finite(tAll$ld10))] <- NA
hist(tAll$ld10)
summary(tAll$ld10,useNA='always')
# smallest tree now had d10 = 1*1.176+1.07 = 2.246 for d10. log10 of this = 0.3514

# We've now changed to using 2018 size data unless missing, and then using 2013 - this comment and 5 lines of code below are from prior analysis. Left here as reminder about how this influenced the U-shaped multinomials. 

# now move diameters forward from 2013 (or 2018, if we use 2013 above), if it's missing in 2018. This can be commented out so that we only use data from one year or the other, and don't mix. I just tried this, for LN.18 analysis - it gets rid of the U-shaped result. So the result is coming from mixing data from plants censuses only in 2013 with those added in 2018. Now we have to figure out why!
# which(is.na(tAll$ld10) & !is.na(tAll$d10.13))
# tAll$ld10[which(is.na(tAll$ld10))] <- log10(tAll$d10.13[which(is.na(tAll$ld10))])
# hist(tAll$ld10)

# create types to switch between types
types <- c('SA','TR')
op=par(mfrow=c(1,2))
for (i in 1:2) {
  ty <- types[i]
  print(hist(tAll$ld10[tAll$Type.18==ty],main=paste('Type',ty)))
}
par(op)

# create squared variable for quadratic analysis
tAll$ld10.2 <- tAll$ld10^2

# check on fates (for all 6945-703 indvs)
table(tAll$Live.18, tAll$fate.18,useNA='always')

# gCrown not properly coded, fixed here using fates
table(tAll$gCrown.18, tAll$fate.18,useNA='always')
tAll$gCrown.18[which(tAll$fate.18 %in% c('LN','LR'))] <- 1
tAll$gCrown.18[which(tAll$fate.18 %in% c('DN','DR'))] <- 0
table(tAll$gCrown.18, tAll$fate.18,useNA='always')

table(tAll$bSprout.18,tAll$fate.18,useNA='always')

# 703 plants have fate=NA in 2018, let's have a look
# these appear to be points added in 2019 and 2020 - that many? - Yes, especially lots of CEACUN!
NA704 <- which(is.na(tAll$fate.18))
tAll[NA704[3],]

# now try saplings and trees separately
op <- par(mfrow=c(1,2))
for (i in 1:2) plot(Live.18~ld10,data=tAll[tAll$Type.18==types[i],],main=paste('Type',types[i]))
par(op)

# or plot together
plot(Live.18~ld10,data=tAll,main=paste('All',types[i]))

# check fate values (for 6241 indvs- excluding the 703 NA's)
table(tAll$Resprout.18,tAll$fate.18)

#### SURVIVAL ANALYSIS - first cut, size only!!
types <- c('TR','SA')

# adjust values here to subset data, for TR and/or SA and fire severity level and species (if needed)
table(tAll$Type.18,useNA='always')
table(tAll$fsLevel,tAll$Plot,useNA='always')

write.csv(tAll,'data/tAll-analyzeData-update.csv')

## XYZABC - run to here 

### MODELING SECTION STARTS HERE
# CREATE tAlls for modeling
tAlls <- tAll[which(tAll$Type.18 %in% types[1:2] & tAll$fsLevel>=0),]
dim(tAll)
dim(tAlls)
table(tAlls$SpCd14)

# NOW DECIDE WHICH SPECIES WE WANT FOR RUNNING MODELS
tAlls <- tAlls[-which(tAlls$SpCd14 %in% c('BACPIL','QUEBER','OTHER')),]
spAs <- spA[-which(spA %in% c('OTHER','BACPIL','QUEBER'))]
dim(tAlls)

# Optional - remove saplings added in 2018, where we might be introducing detection bias towards small survivors
newSap <- which(is.na(tAlls$Year.13) & tAlls$Year.18==2018 & tAlls$Type.18=='SA')
length(newSap)
tAlls$Num[newSap]
table(tAlls$Plot[newSap])
names(tAlls)

# export
expt <- tAlls[newSap,c('Plot','Num','SA.Height_cm.18','SA.BD_cm.18','Basal.Resprout.Height_cm.18','fsLevel','fate.18','Notes.18')]
tail(expt)
write.csv(expt,'data/newSap-export.csv')

# sizes ofnew saplings
summary(tAlls$ld10[newSap])
summary(tAlls$ld10[tAlls$Type.18=='SA'])
summary(tAlls$ld10[tAlls$Type.18=='TR'])

# OPTION: comment this in or out to exercise option
# tAlls <- tAlls[-newSap,]
table(tAlls$SpCd14)
dim(tAlls)

# OPTION: if needed trim data by stem size
#tAlls <- tAlls[which(tAlls$ld10>=c(-0.5)),]
dim(tAlls)

# how many new trees added - includes new plots - not recruitment
newTrees <- which(is.na(tAlls$Year.13) & tAlls$Year.18==2018 & tAlls$Type.18=='TR')
length(newTrees)
table(tAlls$Plot.18[newTrees])

# how many saplings grew into trees
table(tAlls$Type.13,tAlls$Type.18)

# sample size
tAlls$fPlot <- as.factor(tAlls$Plot.18)
nrow(tAlls)
(N <- length(which(!is.na(tAlls$ld10) & !is.na(tAlls$Live.18))))

# START SECTION 'YVAR_MODEL' (search for that name below to see where it ends)
### START MODELING HERE
#set yvalue
# For some outcomes, like NR, need to drop some species
yvalname <- 'gCrown.18'
tAlls$yval <- tAlls[,yvalname]

#fit model
fit0 <- glm(yval~ld10,data=tAlls,family='binomial')
BIC(fit0)

fit1 <- glm(yval~ld10+ld10.2,data=tAlls,family='binomial')
#summary(fit1)
AIC(fit1)
BIC(fit1)

fit1n <- glm(yval~ld10+ld10.2+northness,data=tAlls,family='binomial')
#summary(fit1)
AIC(fit1n)
BIC(fit1n)
coefficients(fit1n)

# we can check all combinations, but for Live.18 and the full data, northness wasn't justified to add to the full model, with FS and species. So for now we're including it, but need to check, for each variable, before reporting in the paper

fit2 <- glmer(yval~ld10+ld10.2+northness+(1|fPlot),data=tAlls,family='binomial')
#summary(fit2)
AIC(fit2)
BIC(fit2)
fitPlots <- rownames(coefficients(fit2)$fPlot)

# now add fire severity
# check whether model can fit fsLevel and random plot factor - converges!
fit3 <- glm(yval~ld10+ld10.2+northness+fsCat,data=tAlls,family='binomial')
#summary(fit3)
AIC(fit3)
BIC(fit3)

fit4 <- glmer(yval~ld10+ld10.2+northness+fsCat+(1|fPlot),data=tAlls,family='binomial')
#summary(fit4)
BIC(fit4)
# YES! (for both types, absp, all fsLevels....)

## MAIN MODEL WE'RE FOCUSING ON
fit5 <- glm(yval~ld10+ld10.2+fsCat+northness+SpCd14,data=tAlls,family='binomial')
BIC(fit5)
coefficients(fit5)

# here's the full model with FS and species, and no northness - check via BIC whether it can be included in final results
fit5x <- glm(yval~ld10+ld10.2+fsCat+SpCd14,data=tAlls,family='binomial')
length(fit5x$residuals)
BIC(fit5x)

#without quadratic?
fit5l <- glm(yval~ld10+fsCat+northness+SpCd14,data=tAlls,family='binomial')
BIC(fit5l)

#fit6 <- glmer(yval~ld10+ld10.2+fsCat+SpCd14+(1|fPlot),data=tAlls,family='binomial')
# DOESN'T CONVERGE combining species and random factor plots

# made newdata for prediction
nvals <- 13
ld10vals <- seq(min(tAlls$ld10,na.rm=T),max(tAlls$ld10,na.rm=T),length.out=nvals)

#DA 1/2/24 - added 100 here, as well as max value 
ld10vals <- c(min(tAlls$ld10,na.rm=T),log10(c(1,1.5,2,5,10,15,20,30,50,75,100)),max(tAlls$ld10,na.rm=T))

ndf <- expand.grid(ld10vals,fitPlots,sort(unique(tAlls$fsCat)),spAs)
names(ndf) <- c('ld10','fPlot','fsCat','SpCd14')
ndf$fPlot <- as.factor(ndf$fPlot)
ndf$fsCat <- as.factor(ndf$fsCat)
ndf$SpCd14 <- as.character(ndf$SpCd14)
ndf$ld10.2 <- ndf$ld10^2
ndf$northness <- 0
ndf$eastness <- 0
head(ndf)
dim(ndf)

nd <- expand.grid(ld10vals,sort(unique(tAlls$fsCat)),spAs)
names(nd) <- c('ld10','fsCat','SpCd14')
nd$SpCd14 <- as.character(nd$SpCd14)
nd$fsCat <- as.factor(nd$fsCat)
nd$ld10.2 <- nd$ld10^2
nd$northness <- 0
head(nd)
dim(nd)

#xx <- data.frame(ld10=seq(min(tAlls$ld10,na.rm=T),max(tAlls$ld10,na.rm=T),length.out=nvals))
#nd2 <- data.frame(fPlot=rep(as.factor(fitPlots),nrow(xx)),ld10=rep(xx$ld10,each=length(fitPlots)))
#nd3 <- data.frame(fPlot=rep(as.factor(fitPlots),nrow(nd2)),ld10=rep(nd2$ld10,nrow(nd2)),fsCat <- )
#rm('xx')


# predict value from fit
nd$predVal0 <- predict(fit0,newdata=nd,type='response')
nd$predVal1 <- predict(fit1,newdata=nd,type='response')
ndf$predVal2 <- predict(fit2,newdata=ndf,type='response')
nd$predVal3 <- predict(fit3,newdata=nd,type='response')
ndf$predVal4 <- predict(fit4,newdata=ndf,type='response')
nd$predVal5 <- predict(fit5,newdata=nd,type='response')
nd$predVal5l <- predict(fit5l,newdata=nd,type='response')
head(nd)

#plot data and predicted values
range(tAlls$ld10,na.rm=T)

plotToFile <- T
getwd()
if (plotToFile) png('figures/outputfig.png',width = 600,height = 400)
plot(tAlls$ld10,tAlls$yval,main=yvalname)
#points(nd$ld10,nd$predVal0,lwd=4) # linear size
#points(nd$ld10,nd$predVal1,lwd=4) # quadratic size
#points(ndf$ld10,nd$predVal2,lwd=4) # qsize + plots only
#points(nd$ld10,nd$predVal3,lwd=4) # qsize + fire levels
#points(ndf$ld10,nd$predVal4,lwd=4) # qsize + fire levels + plots

# comment out next line, and uncomment below if you want to see linear size model instead
points(nd$ld10,nd$predVal5,lwd=4) # qsize + fire levels + species
#points(nd$ld10,nd$predVal5l,lwd=4) # qsize + fire levels + species
if (plotToFile) {
  #system('open figures/outputfig.png')
  dev.off()
}

#rowSel <- which(nd$Species.18=='UMBCAL' & nd$fsCat==1)
#points(nd$ld10[rowSel],nd$predVal5[rowSel],lwd=4,col='red') # qsize + fire levels + species

# visualize mean fire severity - no plot command!
# pAvg <- rep(NA,11)
# i=1
# fs=1
# for (fsv in 1:4) {
#   fs <- fsv-1
#   for (i in 1:nvals) {
#     ld <- ld10vals[i]
#     # comment/uncomment to switch between quadratic (5) and linear (5l)
#     pAvg[i] <- mean(nd$predVal5[which(nd$ld10==ld & nd$fsCat==fs)])
#     #pAvg[i] <- mean(nd$predVal5l[which(nd$ld10==ld & nd$fsCat==fs)])
#   }
#   #print(pAvg)
#   points(ld10vals,pAvg,col='red',lwd=4,type='b')
# }

#Visualize individual species 
par(mfrow=c(1,1))
spname='ARBMEN'

plotSpecies <- function(spname,tmp=tAlls,ndt=nd,ylab=yvalname) {
  tmp <- tmp[which(tmp$SpCd14 %in% spname),]
  plot(tmp$ld10,tmp$yval,xlim=range(tmp$ld10,na.rm=T),main=paste(ylab,spname))
  ndt <- ndt[which(ndt$SpCd14==spname),]
  i=2
  for (i in 1:length(unique(ndt$fsCat))) {
    ndt2 <- ndt[which(ndt$fsCat==(i-1)),]
    lines(ndt2$ld10,ndt2$predVal5)
    #text(ld10.vals[5],ndt$pSurvAll[which(ndt$ld10==ld10.vals[5])],fsLevels.vals[i])
  }
}

#This will plot isoclines of survival at each fire severity for any indv species 
plotSpecies('PSEMEN')
plotSpecies('ARBMEN')
plotSpecies('QUEAGR')
plotSpecies('QUEDOU')
plotSpecies('UMBCAL')
plotSpecies('QUEGAR')
plotSpecies('QUEKEL')
plotSpecies('HETARB')
plotSpecies('ARCMAN')
plotSpecies('AMOCAL')
plotSpecies('FRACAL')

# see predicted values for each species at selected size and fire severity. These can be plotted against bark thickness! Use predVal5 or predVal5l. Others don't have species so they are all the same. ld10=1 is for basal diameter of 10cm. Nice spread among species.
nd[which(nd$ld10==1&nd$fsCat==1),]

### read in bark thickness here and transfer to nd predicted value data.frame
bt <- read.csv('input_data/barkthickness_spmeans.csv')
head(bt)
b2p <- match(nd$SpCd14,bt$names)
nd$t10 <- bt$t10[b2p]
nd$t20 <- bt$t20[b2p]
nd$t30 <- bt$t30[b2p]
nd$t50 <- bt$t50[b2p]

head(nd)
dim(nd)
table(nd$fsCat)
table(10^(nd$ld10))

write.csv(nd,paste('data/',yvalname,'-predValues.csv',sep=''))

# DA 1/2/24 - everything working up to here!!

# explore - NOT MUCH GOING ON!!!
tmp <- nd[which(nd$ld10==log10(30)&nd$fsCat==3),]
plot(predVal5~t30,data=tmp)
text(tmp$t30,tmp$predVal5,labels=tmp$Species.18)

#abline(h=0.5,lty=2) 
#h draws horizontal line at .5, lty -dashed or solid, lwd is line width

## We have an issue with apparent higher survival in small individuals. Let's investigate
#plot(Live.18~ld10,data=tAlls)

# next section not edited to separate TR and SA, not valid!!
if (FALSE) {
  #what is critical basal area to achieve 50% survival?
  (ld50 <- nd$ld10[which(nd$pSurvAll>=0.5)[1]])
  abline(v=ld50,lty=2)
  
  (spRes <- data.frame(species='All',N=N,ld50=ld50,slp=fit$coefficients[2]))
  spRes
  
  
  i <- 3
  for (i in 1:length(spA))
  {
    (species <- spA[i])
    tAlls <- tAll[which(tAll$Species.13==species),]
    N <- length(which(!is.na(tAlls$Basal.Area.13) & !is.na(tAlls$Live.18)))
    fit <- glm(Live.18~ld10,data=tAlls,family='binomial')
    pval <- predict(fit,newdata=nd,type='response')
    nd <- data.frame(nd,pval)
    names(nd)[length(names(nd))] <- paste('pSurv_',species,sep='')
    lines(nd$ld10,nd[,ncol(nd)])
    
    ld50 <- NA
    if (fit$coefficients[2]>0) 
      if (nd[1,ncol(nd)]<0.5) 
        ld50 <- nd$ld10[which(nd[,ncol(nd)]>=0.5)[1]] 
    spRes <- rbind(spRes,c(species,N,ld50,fit$coefficients[2]))
  }
  spRes$ld50 <- as.numeric(spRes$ld50)
  spRes$d50 <- round(10^spRes$ld50,3)
  spRes$ld50 <- round(as.numeric(spRes$ld50),3)
  spRes$slp <- round(as.numeric(spRes$slp),3)
  spRes
  
  plotSP <- function(tAll,species=NULL,xlims=range(tAll$ld10,na.rm=T))
  {
    tAlls <- tAll[which(tAll$Species.13==species),]
    mindbh <- min(tAlls$ld10,na.rm=T)
    maxdbh <- max(tAlls$ld10,na.rm=T)
    plot(tAlls$ld10,tAlls$Live.18,main=species,xlim=xlims)
    fit <- glm(Live.18~ld10,data=tAlls,family='binomial')
    nd <- with(tAll,data.frame(ld10=seq(mindbh,maxdbh,length.out=1001)))
    pval <- predict(fit,newdata=nd,type='response')
    lines(nd$ld10,pval)  
  }
  plotSP(tAll,'ARCMAN')
  plotSP(tAll,'QUEAGR')
  plotSP(tAll,'QUEDOU')
  plotSP(tAll,'ARBMEN')
  plotSP(tAll,'HETARB')
  plotSP(tAll,'PSEMEN')
  plotSP(tAll,'QUEGAR')
  plotSP(tAll,'QUEKEL')
  plotSP(tAll,'UMBCAL')
}

########### Starting with Models########
head(nd)

# 
predSizes <- c(5,10,30,75)
fsl <- 1

# with shrubs
selSize <- 1 #change size 1 above from 10cm to 2cm in order to look at shrubs
barplot(predVal5~SpCd14,data=nd[which(nd$ld10==log10(predSizes[selSize]) & nd$fsCat==fsl),],ylim=c(0,1),main=paste(yvalname,'@ FSLev',fsl,'predval@',predSizes[selSize],'cm dbh'))

# remove shrubs for larger sizes cm plot
# op=par will stack 3 figures (1,3) for horizontal & (3,1) for vertical
op=par(mfrow=c(1,3)) 
for (selSize in 2:4) {
  tmp <- nd[which(!nd$SpCd14 %in% c('HETARB','ARCMAN')),]
  barplot(predVal5~SpCd14,data=tmp[which(tmp$ld10==log10(predSizes[selSize]) & tmp$fsCat==fsl),],ylim=c(0,1),main=paste(yvalname,'@ FSLev',fsl,'predval@',predSizes[selSize],'cm dbh'))
}
par(op)
## END SECTION 'YVAR_MODEL'. THIS IS THE END OF THE ANALYSES AND VISUALIZATIONS THAT START WITH THE SELECTION OF A DEPENDENT VARIABLE.

# 7/7/23 - everything working to here!!!!

# SECTION 'THREE_FATES' - binomial and multinomial compared, for visual purposes
tAlls$Dead.18 <- 1 - tAlls$Live.18
tAlls$fate3.18 <- (-1)
tAlls$fate3.18[which(tAlls$fate.18=='DN')] <- 0
tAlls$fate3.18[which(tAlls$fate.18=='DR')] <- 1
tAlls$fate3.18[which(tAlls$fate.18 %in% c('LN','LR'))] <- 2
table(tAlls$fate3.18) 

(f3Levels <- 0:2)
(f3PlotVals <- c(0.95,1,1.05))
f3Labs <- c('DN','DR','LR+LN')
(f3PlotCols <- c('black','red','green'))

tAlls$f3PlotVals <- f3PlotVals[match(tAlls$fate3.18,f3Levels)]

# comment in or out to select one of these options
FireLevels <- c('Mod+High'); FVals <- 2:3
FireLevels <- c('Low'); FVals <- 1
#FireLevels <- c('None'); FVals <- 0
#FireLevels <- c('Low:High'); FVals <- 1:3 #changes FireSev range to all c(1:3)

#### NEED TO SUBSET BY TYPE HERE
tAllsp <- tAlls[which(tAlls$fsCat %in% FVals & tAlls$Species.18 %in% spA),]

# remove rows with no size data
tAllsp <- tAllsp[which(!is.na(tAllsp$ld10)),]
dim(tAllsp)

# FIT EACH FATE SEPARATELY AND THEN COMPARE MULTINOMIAL
plot(f3PlotVals~ld10,data=tAllsp,col=tAllsp$f3PlotCols,ylim=c(0,1.1),main=paste('Fire Level:',FireLevels))

fit1 <- glm(Dead.18~ld10 +ld10.2,data=tAllsp)
points(fitted(fit1)~tAllsp$ld10,col=f3PlotCols[1])

fit2 <- glm(DR.18~ld10 +ld10.2,data=tAllsp)
points(fitted(fit2)~tAllsp$ld10,col=f3PlotCols[2])

fit3 <- glm(gCrown.18~ld10 +ld10.2,data=tAllsp)
points(fitted(fit3)~tAllsp$ld10,col=f3PlotCols[3])

# NOW FIT MULTINOMIAL
require(nnet)

# MULTINOMIAL - QUADRATIC CAN BE ADDED HERE '+ld10.2' - changes results some
fit.mn <- multinom(as.factor(fate3.18) ~ ld10 + ld10.2, data=tAllsp)
fit.mn 
head(round(fitted(fit.mn),2))
dim(fitted(fit.mn))

plot(f3PlotVals~ld10,data=tAllsp,col=tAllsp$f3PlotCols,ylim=c(0,1.1),main=paste('Fire Level:',FireLevels))
for (i in 1:3) points(tAllsp$ld10,fitted(fit.mn)[,i],col=f3PlotCols[i])

summary(apply(fitted(fit1)[,1:3],1,sum))
#######################

# Full Multinomial model
tAllsp <- tAlls[which(!is.na(tAlls$ld10)),]
dim(tAllsp)

fit.mn <- multinom(as.factor(fate3.18) ~ ld10 + ld10.2 + northness + fsCat + SpCd14, data=tAllsp)
fit.mn
dim(fitted(fit.mn))
for (i in 1:3) plot(fitted(fit.mn)[,i]~tAllsp$ld10,col=f3PlotCols[i],main=paste('Full multinomial',f3Labs[i]))
BIC(fit.mn)
str(fitted(fit.mn))

plot(f3PlotVals~ld10,data=tAllsp,col=tAllsp$f3PlotCols,ylim=c(0,1.1),main='Full multinomial')
selR <- which(tAllsp$Species.13=='QUEAGR' & tAllsp$fsCat==3)
for (i in 1:3) points(tAllsp$ld10,fitted(fit.mn)[,i],col=f3PlotCols[i])
for (i in 1:3) points(tAllsp$ld10[selR],fitted(fit.mn)[selR,i],col= 'purple',pch=19)

# try drop 1
#drop1(fit.mn,c('ld10'))

# drop northness to check if all justified
fit.mnx <- multinom(as.factor(fate3.18) ~ ld10 + ld10.2 + fsCat + Species.18, data=tAllsp)
BIC(fit.mn);BIC(fit.mnx)

plot(f3PlotVals~ld10,data=tAllsp,col=tAllsp$f3PlotCols,ylim=c(0,1.1),main='Full multinomial')
selR <- which(tAllsp$Species.13=='QUEAGR' & tAllsp$fsCat==2)
for (i in 1:3) points(tAllsp$ld10,fitted(fit.mnx)[,i],col=f3PlotCols[i])
for (i in 1:3) points(tAllsp$ld10[selR],fitted(fit.mnx)[selR,i],col= 'purple',pch=19)
# northness supported

# END CLEAN CODE HERE!!! redundant stuff happening below


###old figures- this doesn't work anymore###

# spA
# png('figures/logit1.png',width = 800, height = 1200)
# op=par(mfrow=c(4,2))
# plotSpecies('AMOCAL')
# plotSpecies('ARBMEN')
# plotSpecies('ARCMAN')
# plotSpecies('FRACAL')
# par(op)
# dev.off()
# 
# png('figures/logit2.png',width = 800, height = 1200)
# op=par(mfrow=c(4,2))
# plotSpecies('HETARB')
# plotSpecies('PSEMEN')
# plotSpecies('UMBCAL')
# par(op)
# dev.off()
# 
# png('figures/logit3.png',width = 800, height = 1200)
# op=par(mfrow=c(4,2))
# plotSpecies('QUEAGR')
# plotSpecies('QUEDOU')
# plotSpecies('QUEGAR')
# plotSpecies('QUEKEL')
# par(op)
# dev.off()

### full model with random plot
#fit <- glmer(Live.18~ld10 + FireSev + Species.13 + (1 |Plot.13),data=tAlls,family='binomial')

######## TOPKILL ANALYSIS
# current scoring has dead trees as topkilled. Change so topkill is only for those that are alive
tAll$TopkillLive.18 <- 0
tAll$TopkillLive.18[which(tAll$Topkill.18==1 & tAll$Live.18==1)] <- 1

table(tAll$Dead.18,tAll$TopkillLive.18,tAll$gCrown.18)
tAll[which(tAll$gCrown.18==0 & tAll$TopkillLive.18==0 & tAll$Dead.18==0),]
# TWO INDIVIDUALS WITH PROBLEM DATA: 3415, 4437 (already identified those above - if they've been fixed and don't show up at this point, delete this line)

tAlla <- tAll[which(tAll$FireSev>100),]

dim(tAlla)
names(tAlla)
head(tAlla)
tAlla$ld10.2 <- tAlla$ld10^2

op=par(mfrow=c(1,1))
plot(Dead.18~ld10,data=tAlla)
fitD <- glm(Dead.18~ld10+ld10.2,data=tAlla,family='binomial')
fitD
nd <- data.frame(ld10=seq(min(tAlla$ld10,na.rm=T),max(tAlla$ld10,na.rm=T),length.out=101))
nd$ld10.2 <- nd$ld10^2
head(nd)
nd$pDead <- predict(fitD,nd,type='response')
lines(pDead~ld10,data=nd,lwd=2,col='black')

#plot(Topkill.18~ld10,data=tAlla)
fitT <- glm(TopkillLive.18~ld10+ld10.2,data=tAlla,family='binomial')
fitT
nd$pTopKill <- predict(fitT,nd,type='response')
lines(pTopKill~ld10,data=nd,lwd=2,col='red')

#plot(gCrown.18~ld10,data=tAlla)
fitG2 <- glm(gCrown.18~ld10+ld10.2,data=tAlla,family='binomial')
fitG2
BIC(fitG2)
fitG1 <- glm(gCrown.18~ld10,data=tAlla,family='binomial')
fitG1
BIC(fitG1)

nd$pGreen <- predict(fitG1,nd,type='response')
lines(pGreen~ld10,data=nd,lwd=2,col='green')

par(op)

#### NOW TRY MULTINOMIAL### SKIP THIS SECTION-###
require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(tAll[,c('ld10','PFstatus','Species.13','FireSev')],1,allNotNA)
table(rComp)
tAlla <- tAll[rComp,]

fit1 <- multinom(PFstatus ~ ld10,data=tAlla)
fit1
head(round(fitted(fit1),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit1)[,1],col='red')
points(tAlla$ld10,fitted(fit1)[,2],col='gray')
points(tAlla$ld10,fitted(fit1)[,3],col='green')


##########
head(tAll)
dim(tAll)
tAll$dupStatusCheck <- apply(tAll[,c('Dead.18','TopkillLive.18','gCrown.18')],1,sum,na.rm=T)
table(tAll$dupStatusCheck)
## same 2 as above

###end skip section ###

#### PROVISIONALLY ASSIGN TO THREE CLASSES
tAll$PFstatus <- (-1)
tAll$PFstatus[which(tAll$Live.18==0)] <- 0
tAll$PFstatus[which(tAll$TopkillLive.18==1)] <- 1
tAll$PFstatus[which(tAll$Topkill.18==0 & tAll$gCrown.18==1)] <- 2
table(tAll$PFstatus)

# now assign all remaing NAs to dead - TEMP STEP
tAll$PFstatus[which(tAll$PFstatus==(-1))] <- 0
table(tAll$PFstatus)

## NOW TRY MULTINOMIAL
# setup discrete FireSev
summary(tAll$FireSev)
tAll$dFS <- cut(tAll$FireSev,breaks = c(-200,35,130,298,1000))
table(tAll$dFS)

require(nnet)
allNotNA <- function(x) all(!is.na(x))

# Size only
rComp <- apply(tAll[,c('ld10','PFstatus','Species.13','FireSev')],1,allNotNA)
table(rComp)
tAlla <- tAll[rComp,]

fit1 <- multinom(PFstatus ~ ld10,data=tAlla)
fit1
head(round(fitted(fit1),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit1)[,1],col='red')
points(tAlla$ld10,fitted(fit1)[,2],col='gray')
points(tAlla$ld10,fitted(fit1)[,3],col='green')

#Size and Fire Sev 
fit2 <- multinom(PFstatus ~ ld10 + dFS,data=tAlla)
fit2
BIC(fit2)
head(round(fitted(fit2),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit2)[,1],col='red')
points(tAlla$ld10,fitted(fit2)[,2],col='grey')
points(tAlla$ld10,fitted(fit2)[,3],col='green')

# Add Species
fit3 <- multinom(PFstatus ~ ld10 + dFS + Species.13,data=tAlla)
fit3
BIC(fit3)
head(round(fitted(fit3),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit3)[,1],col='red')
points(tAlla$ld10,fitted(fit3)[,2],col='grey')
points(tAlla$ld10,fitted(fit3)[,3],col='green')

# Add Species * size interaction
fit3x <- multinom(PFstatus ~ ld10 + dFS + Species.13 + dFS:ld10,data=tAlla)
fit3x
BIC(fit3x)
head(round(fitted(fit3x),2))

plot(tAlla$ld10,tAlla$Live.18)
points(tAlla$ld10,fitted(fit3x)[,1],col='red')
points(tAlla$ld10,fitted(fit3x)[,2],col='black')
points(tAlla$ld10,fitted(fit3x)[,3],col='green')

# COMPARE BIC
BIC(fit1)
BIC(fit2)
BIC(fit3)
BIC(fit3x)

# plot for selected values
rsel <- which(tAlla$Species.13=='PSEMEN' & as.numeric(tAlla$dFS)>0)
plot(tAlla$ld10[rsel],tAlla$PFstatus[rsel])
points(tAlla$ld10[rsel],fitted(fit3x)[rsel,1],col='red',pch=19)
points(tAlla$ld10[rsel],fitted(fit3x)[rsel,2],col='black',pch=19)
points(tAlla$ld10[rsel],fitted(fit3x)[rsel,3],col='green',pch=19)

## logit model and plot
d=tAll
xcn='ld10'
ycn='Live.18'
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
logit2Plot(tAll,'ld10','Live.18')




#### MOVED DOWN FROM JUST BEFORE MULTINOMIAL
if (FALSE) {
  dim(tAllsp)
  
  # code to run a binomial model and plot response curve with data
  # MORTALITY
  nd <- with(tAllsp,data.frame(ld10=seq(min(tAllsp$ld10,na.rm=T),max(tAllsp$ld10,na.rm=T),length.out=101)))
  nd$ld10.2 <- nd$ld10^2
  
  #plot(tAllsp$ld10,tAllsp$Dead.18,xlim=tAllld10.range)
  fit1 <- glm(Dead.18~ld10,data=tAllsp,family='binomial')
  fit2 <- glm(Dead.18~ld10+ld10.2,data=tAllsp,family='binomial')
  BIC(fit1)
  BIC(fit2)
  nd$pMortality <- predict(fit1,newdata=nd,type='response')
  #lines(nd$ld10,nd$pMortality)
  
  #TOPKILL WITH RESPROUT
  #plot(tAllsp$ld10,tAllsp$TB,xlim=tAllld10.range)
  fit <- glm(DR.18~ld10+ld10.2 ,data=tAllsp,family='binomial')
  summary(fit)
  nd$pDR <- predict(fit,newdata=nd,type='response')
  #lines(nd$ld10,nd$pTB)
  
  #GREEN CROWN
  #plot(tAllsp$ld10,tAllsp$gCrown.18,xlim=tAllld10.range)
  fit1 <- glm(gCrown.18~ld10,data=tAllsp,family='binomial')
  fit2 <- glm(gCrown.18~ld10+ld10.2,data=tAllsp,family='binomial')
  summary(fit1)
  nd$pGCrown <- predict(fit1,newdata=nd,type='response')
  #lines(nd$ld10,nd$pGCrown)
  
  #plot all three (change main from selSpecies to ".." to alter main title)
  plot(tAllsp$ld10,tAllsp$PFsPlotVals,col=tAllsp$PFsPlotCols,pch=19,ylim=c(-0.05,1.05),xlim=tAllld10.range,main=paste("ALL_SP",FireLevels))
  points(tAllsp$ld10,rep(-0.05,length(tAllsp$ld10)))
  lines(nd$ld10,nd$pGCrown,col='green')
  lines(nd$ld10,nd$pDR,col='red')
  lines(nd$ld10,nd$pMortality)
  
  # does the sum of the three binomials for these three exclusive fates sum to 1?
  nd$pTOT <- apply(nd[,c('pGCrown','pDR','pMortality')],1,sum)
  summary(nd$pTOT)
  #reset plotting window with one panel
  par(mfrow=c(1,1))
  
}

