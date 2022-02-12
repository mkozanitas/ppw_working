# RUN prepareData and examineData again if csv's have been changed, or get.indv.data() has been updated in PW_functions_local.R
rm(list=ls())

lba2d <- function(x) 2*sqrt((10^x)/pi)
ba2d <- function(x) 2*sqrt((x)/pi)
d2ba <- function(x) pi*(x/2)^2
d2lba <- function(x) log10(pi*(x/2)^2)

all.id <- readRDS('data/allid-nodups.Rdata')
str(all.id)

spNames <- read.csv('data/all-spp-names.csv')
head(spNames)
allIndv <- readRDS('data/allIndv.Rdata')
head(allIndv)

## get percent survival by species and type, from 13-18
fst <- data.frame(Species=rep(spNames$x,each=2),Type=rep(c('SA','TR'),nrow(spNames)),N13=NA,N18=NA,S13.18=NA)
head(fst)                  

# select which individuals to eliminate - this can be done separately for each analysis
# take any individuals where plot and species match for the first two years
nrow(allIndv)
plot.ok <- allIndv$Num[which(allIndv$P13==allIndv$P18)] 
spec.ok <- allIndv$Num[which(allIndv$S13==allIndv$S18)] 
ps.ok <- intersect(plot.ok,spec.ok)
length(ps.ok)

# not missed in 2018 survey
nm <- all.id[[1]]$Num[which(all.id[[1]]$NA18==0)]
head(nm)
length(nm)

t1 <- all.id[[1]][which(intersect(ps.ok,nm) %in% all.id[[1]]$Num),]
t1 <- t1[order(t1$Num),]
nrow(t1)
t2 <- all.id[[2]][which(intersect(ps.ok,nm) %in% all.id[[2]]$Num),]
nrow(t2)
t2 <- t2[order(t1$Num),]

# UNMMERGED - ANALYSES BASED ON SHARED NUMS
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
head(fst)
fst$S13.18[which(fst$N13!=0)] <- fst$N18[which(fst$N13!=0)]/fst$N13[which(fst$N13!=0)]
fst[order(fst$N13,decreasing = T),]

sum(fst$N13)
sum(fst$N18)
sum(fst$N18)/sum(fst$N13)

fst
barplot(t(as.matrix(fst[,c('N13','N18')])))

# NOW MERGE FOR ANALYSES BY INDIV
t12 <- merge(t1,t2,'Num',all = T)
dim(t12)
names(t12)
head(which(is.na(t12$Plot.x)))

# CALC DBH - easier to connect with our field data and knowledge
t12$dbh <- ba2d(t12$Basal.Area.x)
summary(t12$dbh)
plot(log10(t12$dbh),t12$Live.y)

## USE LOG DBH for analysis
t12$ldbh <- log10(t12$dbh)
hist(t12$ldbh)

# **** mixing SA and TR - need to work on diameter conversion to get this right ***
N <- length(which(!is.na(t12$ldbh) & !is.na(t12$Live.y)))
fit <- glm(Live.y~ldbh,data=t12,family='binomial')
summary(fit)

nd <- with(t12,data.frame(ldbh=seq(-0.5,2,length.out=1001)))
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)

plot(t12$ldbh,t12$Live.y)
lines(nd$ldbh,nd$pSurvAll,lwd=4)
abline(h=0.5,lty=2)

#what is critical basal area to achieve 50% survival?
ld50 <- nd$ldbh[which(nd$pSurvAll>=0.5)[1]]
abline(v=ld50,lty=2)
(spRes <- data.frame(species='All',N=N,d50=10^ld50,slp=fit$coefficients[2]))

## now run by species for species with lots of data
spN <- table(t12$Species.x)
spN[order(spN)]
(spA <- names(spN)[which(spN>=25)])
spA <- spA[-which(spA=='QUEBER')] # drop QUEBER - none died!

i <- 3
for (i in 1:length(spA))
{
  (species <- spA[i])
  tmp <- t12[which(t12$Species.x==species),]
  N <- length(which(!is.na(tmp$Basal.Area.x) & !is.na(tmp$Live.y)))
  fit <- glm(Live.y~ldbh,data=tmp,family='binomial')
  pval <- predict(fit,newdata=nd,type='response')
  nd <- data.frame(nd,pval)
  names(nd)[length(names(nd))] <- paste('pSurv_',species,sep='')
  lines(nd$ldbh,nd[,ncol(nd)])
  
  ld50 <- NA
  if (fit$coefficients[2]>0) 
    if (nd[1,ncol(nd)]<0.5) 
      ld50 <- nd$ldbh[which(nd[,ncol(nd)]>=0.5)[1]] 
  spRes <- rbind(spRes,c(species,N,10^ld50,fit$coefficients[2]))
}
spRes$d50 <- round(as.numeric(spRes$d50),3)
spRes$slp <- round(as.numeric(spRes$slp),3)
spRes

plotSP <- function(t12,species=NULL)
{
  tmp <- t12[which(t12$Species.x==species),]
  plot(tmp$ldbh,tmp$Live.y)
  fit <- glm(Live.y~ldbh,data=tmp,family='binomial')
  nd <- with(t12,data.frame(ldbh=seq(-0.5,2,length.out=1001)))
  pval <- predict(fit,newdata=nd,type='response')
  lines(nd$ldbh,pval)  
}
plotSP(t12,'QUEGAR')


## Full model with species
t12s <- t12[which(t12$Species.x %in% spA),]
(N <- length(which(!is.na(t12s$ldbh) & !is.na(t12s$Live.y))))

fit <- glm(Live.y~ldbh * Species.x,data=t12s,family='binomial')
summary(fit)

nd <- with(t12s,data.frame(ldbh=seq(-0.5,2,length.out=101)))
nd <- with(t12s,data.frame(ldbh=0.5,Species.x=unique(Species.x)))
nd$pSurvAll <- predict(fit,newdata=nd,type='response')
head(nd)
barplot(nd$pSurvAll~nd$Species.x)
