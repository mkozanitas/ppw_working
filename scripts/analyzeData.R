# RUN prepareData and examineData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R
rm(list=ls())

all.id <- readRDS('data/allid-nodups.Rdata')
str(all.id)

spNames <- read.csv('data/all-spp-names.csv')
head(spNames)
allIndv <- readRDS('data/allIndv.Rdata')

## get percent survival by species and type, from 13-18
fst <- data.frame(Species=rep(spNames$x,each=2),Type=rep(c('SA','TR'),nrow(spNames)),N13=NA,N18=NA,S13.18=NA)
head(fst)                  

# select which individuals to eliminate - this can be done separately for each analysis
# take any individuals where plot and species match for the first two years
plot.ok <- allIndv$Num[which(allIndv$P13==allIndv$P18)] 
spec.ok <- allIndv$Num[which(allIndv$S13==allIndv$S18)] 
ps.ok <- intersect(plot.ok,spec.ok)
length(ps.ok)

# not missed in 2018 survey
nm <- all.id[[1]]$Num[which(all.id[[1]]$NA18==0)]
length(nm)

t1 <- all.id[[1]][intersect(ps.ok,nm),]
t2 <- all.id[[2]][intersect(ps.ok,nm),]

i=1
for (i in 1:nrow(fst))
{
  sp <- fst$Species[i]
  ty <- fst$Type[i]
  init <- t1$Num[which(t1$Species==sp & t1$Type==ty & t1$Live==1)]
  fst$N13[i] <- length(init)
  fin <- intersect(init,t2$Num[which(t2$Live==1)])
  # fin is the number of original ty surviving, whether or not they transitioned from SA->TR (or the other way!)
  fst$N18[i] <- length(fin)
}
fst$S13.18[which(fst$N13!=0)] <- fst$N18[which(fst$N13!=0)]/fst$N13[which(fst$N13!=0)]
fst[order(fst$N13,decreasing = T),]

sum(fst$N13)
sum(fst$N18)
sum(fst$N18)/sum(fst$N13)

