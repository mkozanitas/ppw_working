options(scipen=999)
library(brms)

d <- read.csv("UMBCAL.csv")

spSel='UMBCAL'
spName <- spSel
table(d$Species)
fs='low-medium' #'all','low-medium','drop-high'
logt=T

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
d <- d[complete.cases(d$fsCat2,d$d10.17,d$yvar),]
table(d$Plot)
dim(d)

d$fatefac <- factor(d$fate3.18)
table(d$fatefac)

multifit1 <- brm(fatefac ~ s(d10.17, by=fsCat) + (1|Plot) + (1|iNum.13), data=d,
                 family="categorical", chains = 2, cores = 2, seed=726, 
                 #backend="cmdstanr"
                 )
summary(multifit1)
conditional_smooths(multifit1)
conditional_effects(multifit1, categorical=TRUE)

