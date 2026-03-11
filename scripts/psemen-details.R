## the doug fir problem
rm(list=ls())
source('scripts/31.functionsForAnalysis.R')
spName <- 'PSEMEN'
mod <- 'brm'
fit.type <- c('MN.Quad','MN.Splk3','MN.Splk6','MN.Splk20')
dd <- readRDS(paste(local.dir,'/',paste(mod,spName,'dd.rds',sep='.'),sep=''))
dim(dd)

table(dd$Survey,dd$Type.17.4)
range(dd$d10.17,na.rm=T)

hb <- seq(-0.6,2.4,by=0.2)
pt <- hist(dd$d10.17[which(dd$Survey=='Plot')],breaks=hb)
ht <- hist(dd$d10.17[which(dd$Survey=='Hect')],breaks=hb)

ptl <- hist(dd$d10.17[which(dd$Survey=='Plot' & dd$Live.18==1)],breaks=hb)
htl <- hist(dd$d10.17[which(dd$Survey=='Hect' & dd$Live.18==1)],breaks=hb)

pd <- data.frame(pt$counts,ht$counts,ptl$counts,htl$counts)
head(pd)
pd$pt.psurv <- pd$ptl.counts/pd$pt.counts
pd$ht.psurv <- pd$htl.counts/pd$ht.counts
pd$psurv <- (pd$ptl.counts+pd$htl.counts)/(pd$pt.counts + pd$ht.counts)
pd

# for low/med
table(dd$fac.fsCat)
ddlm <- dd[which(dd$fac.fsCat=='12.LM'),]
pt <- hist(ddlm$d10.17[which(ddlm$Survey=='Plot')],breaks=hb)
ht <- hist(ddlm$d10.17[which(ddlm$Survey=='Hect')],breaks=hb)

ptl <- hist(ddlm$d10.17[which(ddlm$Survey=='Plot' & ddlm$Live.18==1)],breaks=hb)
htl <- hist(ddlm$d10.17[which(ddlm$Survey=='Hect' & ddlm$Live.18==1)],breaks=hb)

pdlm <- data.frame(pt$counts,ht$counts,ptl$counts,htl$counts)
head(pdlm)
pdlm$pt.psurv <- pdlm$ptl.counts/pdlm$pt.counts
pdlm$ht.psurv <- pdlm$htl.counts/pdlm$ht.counts
pdlm$psurv <- (pdlm$ptl.counts+pdlm$htl.counts)/(pdlm$pt.counts + pdlm$ht.counts)
pdlm

# for high
table(dd$fac.fsCat)
ddh <- dd[which(dd$fac.fsCat=='3.H'),]
pt <- hist(ddh$d10.17[which(ddh$Survey=='Plot')],breaks=hb)
ht <- hist(ddh$d10.17[which(ddh$Survey=='Hect')],breaks=hb)

ptl <- hist(ddh$d10.17[which(ddh$Survey=='Plot' & ddh$Live.18==1)],breaks=hb)
htl <- hist(ddh$d10.17[which(ddh$Survey=='Hect' & ddh$Live.18==1)],breaks=hb)

pdh <- data.frame(pt$counts,ht$counts,ptl$counts,htl$counts)
pdh$pt.psurv <- pdh$ptl.counts/pdh$pt.counts
pdh$ht.psurv <- pdh$htl.counts/pdh$ht.counts
pdh$psurv <- (pdh$ptl.counts+pdh$htl.counts)/(pdh$pt.counts + pdh$ht.counts)
pdh

# for unburned
table(dd$fac.fsCat)
ddu <- dd[which(dd$fac.fsCat=='0.U'),]
pt <- hist(ddu$d10.17[which(ddu$Survey=='Plot')],breaks=hb)
ht <- hist(ddu$d10.17[which(ddu$Survey=='Hect')],breaks=hb)

ptl <- hist(ddu$d10.17[which(ddu$Survey=='Plot' & ddu$Live.18==1)],breaks=hb)
htl <- hist(ddu$d10.17[which(ddu$Survey=='Hect' & ddu$Live.18==1)],breaks=hb)

pdu <- data.frame(pt$counts,ht$counts,ptl$counts,htl$counts)
pdu$pt.psurv <- pdu$ptl.counts/pdu$pt.counts
pdu$ht.psurv <- pdu$htl.counts/pdu$ht.counts
pdu$psurv <- (pdu$ptl.counts+pdu$htl.counts)/(pdu$pt.counts + pdu$ht.counts)
pdu

plot(pdu$psurv,ylim=c(0,1),type='b',pch=19)
points(pdlm$psurv,type='b',pch=19,col='orange')
points(pdh$psurv,type='b',pch=19,col='red')

plot(pdu$pt.psurv,ylim=c(0,1),type='b',pch=19)
points(pdlm$pt.psurv,type='b',pch=19,col='orange')
points(pdh$pt.psurv,type='b',pch=19,col='red')

plot(pd$pt.counts+pd$ht.counts,type='b',pch=19)
points(pd$pt.counts,type='b',pch=19,col='blue')

table(dd$Plot,dd$Survey)
