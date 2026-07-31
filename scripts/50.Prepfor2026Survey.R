## This is final analysis script with code in order required to reproduce results in Results section of manuscript

rm(list=ls())
source('scripts/31.functionsForAnalysis.R')
od <- '/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood Drive/CalFire2026_resurvey'

#require(MuMln)

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

# Are all species in spAtt?
spAtt <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/SpeciesData/all-spp-attributes.csv")
(tAllspec <- sort(unique(tAll$Species)))
tAllspec[which(!tAllspec %in% spAtt$OrigSpecies)]
length(which(tAll$Species=='QUEDEC'))

tAll$FuncGroup <- spAtt$FuncGroup[match(tAll$Species,spAtt$OrigSpecies)]
tAll$FuncGroup[which(tAll$Species=='QUEDEC')] <- 'WHTO'
table(tAll$FuncGroup,useNA='a')

tAll$SampGroup <- tAll$FuncGroup
tAll$SampGroup[which(tAll$Species=='HETARB')] <- 'RShrub.Samp'

pAll <- tAll[which(tAll$Survey=='Plot'),]

# Archive so you can come back to this point
pAll.Archive <- pAll

# Reset if needed
pAll <- pAll.Archive

# Best estimate of status as of 2020
sink(paste(od,'PlotDataSummary.txt',sep='/'),F,split=T)
#print('var names');names(pAll)
print('dim pAll');dim(pAll)

print("table(pAll$Live.17,pAll$fate3.18,useNA='a')")
table(pAll$Live.17,pAll$fate3.18,useNA='a')

print("table(pAll$fate3.18,pAll$fate3.19,useNA='a')")
table(pAll$fate3.18,pAll$fate3.19,useNA='a')

print("table(pAll$fate3.19,pAll$fate3.20,useNA='a')")
table(pAll$fate3.19,pAll$fate3.20,useNA='a')

# which plots were surveyed in 2020
S20plots <- names(table(pAll$Plot[which(!is.na(pAll$fate3.20))]))

# first survey
pAll$FirstSurvey <- NA
pAll$FirstSurvey[which(pAll$Live.13==1)] <- 'F13'
pAll$FirstSurvey[which(is.na(pAll$FirstSurvey) & !is.na(pAll$fate3.18))] <- 'F18'
pAll$FirstSurvey[which(is.na(pAll$FirstSurvey) & !is.na(pAll$fate3.19))] <- 'F19'
pAll$FirstSurvey[which(is.na(pAll$FirstSurvey) & !is.na(pAll$fate3.20))] <- 'F20'

print("table(pAll$FirstSurvey,useNA='a')")
table(pAll$FirstSurvey,useNA='a')

# last year when fate was not NA
print('L19.R refers to plants that were last seen in 2019 and the plot WAS surveyed in 2020')
print('L19.S refers to plants that were last seen in 2019 and the plot WAS NOT surveyed in 2020')
pAll$LastSurvey <- NA
pAll$LastSurvey[which(!is.na(pAll$fate3.20))] <- 'L20'
pAll$LastSurvey[which(is.na(pAll$LastSurvey) & !is.na(pAll$fate3.19) & (pAll$Plot %in% S20plots))] <- 'L19.R' #19.R means last seen in 2019, and plot WAS surveyed in 2020
pAll$LastSurvey[which(is.na(pAll$LastSurvey) & !is.na(pAll$fate3.19) & (!pAll$Plot %in% S20plots))] <- 'L19.S' #19.S means last seen in 2019, and plot WAS NOT surveyed in 2020
pAll$LastSurvey[which(is.na(pAll$LastSurvey) & !is.na(pAll$fate3.18))] <- 'L18'

print("table(pAll$LastSurvey,useNA='a')")
table(pAll$LastSurvey,useNA='a')

print("table(pAll$FirstSurvey,pAll$LastSurvey,useNA='a')")
table(pAll$FirstSurvey,pAll$LastSurvey,useNA='a')

#table(pAll$Species,pAll$LastSurvey)

pAll$SampGroup <- 'Z'
pAll$SampGroup[which(pAll$FuncGroup %in% c('EHRO','NS.Con','WHTO'))] <- 'A'
pAll$SampGroup[which(pAll$Species=='HETARB')] <- 'A'
pAll$SampGroup[which(pAll$Species=='ARCMAN')] <- 'A'
pAll$SampGroup[which(pAll$Species=='TORCAL')] <- 'Z'

print("table(pAll$Func,pAll$LastSurvey,pAll$SampGroup)")
table(pAll$Func,pAll$LastSurvey,pAll$SampGroup)

print("Number to sample by plot")
table(pAll$Plot[which(pAll$SampGroup=='A')],pAll$FuncGroup[which(pAll$SampGroup=='A')])

# How many total plants to sample
print("total number to sample in 2026")
length(which(pAll$LastSurvey %in% c('L19.S','L20') & (pAll$SampGroup=='A')))
table(pAll$LastSurvey,pAll$SampGroup)

sink()

# trial data sheets for 3 plots
# DBH (tree) or d10 (sap) when last seen
# fate when last seen

names(pAll) 
pAll$LastType <- NA
pAll$LastDiam <- NA
pAll$LastSapHt <- NA
pAll$LastRspN <- NA
pAll$LastRspHt <- NA
pAll$LastFate4 <- NA
pAll$Quad <- NA
pAll$X <- NA
pAll$Y <- NA

DBHcols <- grep('DBH',names(pAll))
SAHcols <- grep('SA.Height',names(pAll))
fate4cols <- grep('fate4',names(pAll))

pAll$TP <- 0
pAll$TP[which(pAll$Tag.Pulled.18!=0 | pAll$Tag.Pulled.19!=0 | pAll$Tag.Pulled.19!=0)] <- 1
table(pAll$TP)

r20 <- which(pAll$LastSurvey=='L20')
r19 <- which(pAll$LastSurvey %in% c('L19.R','L19.S'))
r18 <- which(pAll$LastSurvey=='L18')
pAll[r20,c('LastType','LastDiam','LastSapHt','LastRspHt','LastRspN','LastFate4')] <- pAll[r20,c('Type.20','DBH_cm.20','SA.Height_cm.20','Basal.Resprout.Height_cm.20','Basal.Resprout.Count.20','fate4.20')]
pAll[r19,c('LastType','LastDiam','LastSapHt','LastRspHt','LastRspN','LastFate4')] <- pAll[r19,c('Type.19','DBH_cm.19','SA.Height_cm.19','Basal.Resprout.Height_cm.19','Basal.Resprout.Count.19','fate4.19')]
pAll[r18,c('LastType','LastDiam','LastSapHt','LastRspHt','LastRspN','LastFate4')] <- pAll[r18,c('Type.18','DBH_cm.18','SA.Height_cm.18','Basal.Resprout.Height_cm.18','Basal.Resprout.Count.18','fate4.18')]

length(which(is.na(pAll$LastDiam)))
srows <- which(is.na(pAll$LastDiam) & !is.na(pAll$DBH_cm.19))
pAll$LastDiam[srows] <- pAll$DBH_cm.19[srows]
srows <- which(is.na(pAll$LastDiam) & !is.na(pAll$DBH_cm.18))
pAll$LastDiam[srows] <- pAll$DBH_cm.18[srows]
srows <- which(is.na(pAll$LastDiam) & !is.na(pAll$DBH_cm.17))
pAll$LastDiam[srows] <- pAll$DBH_cm.17[srows]

length(which(is.na(pAll$LastSapHt)))
srows <- which(is.na(pAll$LastDiam) & !is.na(pAll$SA.Height_cm.19))
pAll$LastDiam[srows] <- pAll$SA.Height_cm.19[srows]
srows <- which(is.na(pAll$LastDiam) & !is.na(pAll$SA.Height_cm.18))
pAll$LastDiam[srows] <- pAll$SA.Height_cm.18[srows]
srows <- which(is.na(pAll$LastDiam) & !is.na(pAll$SA.Height_cm.17))
pAll$LastDiam[srows] <- pAll$SA.Height_cm.17[srows]
length(which(is.na(pAll$LastSapHt)))
length(which(is.na(pAll$LastDiam) & is.na(pAll$LastSapHt)))

pAll$LastLive <- 'L'
pAll$LastLive[which(pAll$LastFate=='DN')] <- 'X'
table(pAll$LastLive)

r20 <- which(pAll$FirstSurvey=='F20')
r19 <- which(pAll$FirstSurvey=='F19')
r18 <- which(pAll$FirstSurvey=='F18')
r13 <- which(pAll$FirstSurvey=='F13')
pAll[r13,c('Quad','X','Y')] <- pAll[r13,c('Quad.13','X_cm.13','Y_cm.13')]
pAll[r18,c('Quad','X','Y')] <- pAll[r18,c('Quad.18','X_cm.18','Y_cm.18')]
pAll[r19,c('Quad','X','Y')] <- pAll[r19,c('Quad.19','X_cm.19','Y_cm.19')]
pAll[r20,c('Quad','X','Y')] <- pAll[r20,c('Quad.20','X_cm.20','Y_cm.20')]

qcols <- grep('Quad.',names(pAll))
xcols <- grep('X_cm.',names(pAll))
ycols <- grep('Y_cm.',names(pAll))

Qlist <- c('A1','A2','A3','A4','B1','B2','B3','B4','C1','C2','C3','C4','D1','D2','D3','D4')
which(!pAll$Quad %in% Qlist)

i=3
for (i in 2:length(qcols))
{
  dq <- which(is.na(pAll[,qcols[i-1]]) & !is.na(pAll[,qcols[i]]) | pAll[,qcols[i]] != pAll[,qcols[i-1]])
  # in 2020 several plots were entered as xI rather than x1 - where x is a column letter
  # in 2019 one plot was entered as 3B rather than B3. Once these are excluded, this loop seems to fix all other quad changes
  dq <- intersect(dq,which(pAll[,qcols[i]] %in% Qlist))
  print(dq)
  pAll[dq,'Quad'] <- pAll[dq,qcols[i]]
  
  xq <- which(is.na(pAll[,xcols[i-1]]) & !is.na(pAll[,xcols[i]]) | pAll[,xcols[i]] != pAll[,xcols[i-1]])
  pAll$X[xq] <- pAll[xq,xcols[i]]
  
  yq <- which(is.na(pAll[,ycols[i-1]]) & !is.na(pAll[,ycols[i]]) | pAll[,ycols[i]] != pAll[,ycols[i-1]])
  pAll$Y[yq] <- pAll[yq,ycols[i]]
}

badQ <- which(!pAll$Quad %in% c('A1','A2','A3','A4','B1','B2','B3','B4','C1','C2','C3','C4','D1','D2','D3','D4'))
pAll[badQ,c(qcols,which(names(pAll)%in%c('Quad','TP','fate4.18')))]
badX <- which(pAll$X>500)
badY <- which(pAll$Y>500)
pAll[union(badX,badY),c(which(names(pAll)%in%c('Plot','Num','X','Y')))]


NumN <- table(pAll$Num)
print(NumN[which(NumN>1)])
pAll[which(pAll$Num==3481),]

# Make tables used by field crew for resampling
sRows <- which(pAll$SampGroup=='A' & pAll$TP==0)
sRows <- which(pAll$TP==0)

length(sRows)
Sdat <- pAll[sRows,c('SampGroup','Plot','Num','Species','Quad','X','Y','FirstSurvey','LastSurvey','LastType','LastDiam','LastSapHt','LastRspN','LastRspHt','LastFate4','LastLive')]

#reorder
fate.order <- data.frame(fate4=c('DN','DR','LN','LR'),order=c('B','A','A','A'))
# DL = sort dead last
Sdat$DL <- fate.order$order[match(Sdat$LastFate4,fate.order$fate4)]
Sdat2 <- Sdat[order(Sdat$Plot,Sdat$SampGroup,Sdat$Quad,Sdat$DL,Sdat$Num),]
head(Sdat2)

table(Sdat2$Plot,Sdat2$SampGroup)

write.csv(Sdat2,paste(od,'AllPlots2026.csv',sep='/'))

pRows <- which(Sdat2$Plot=='PPW1340')
write.csv(Sdat2[pRows,],paste(od,'PPW1340.csv',sep='/'))


pRows <- which(Sdat$Plot=='PPW1341')
write.csv(Sdat[pRows,],paste(od,'PPW1341.csv',sep='/'))


pRows <- which(Sdat$Plot=='PPW1335')
write.csv(Sdat[pRows,],paste(od,'PPW1335.csv',sep='/'))

pRows <- which(Sdat$Plot=='PPW1345')
write.csv(Sdat[pRows,],paste(od,'PPW1345.csv',sep='/'))

# Data for Lauren
# Tree status
# Species
# Live/Dead
# DBH
# Resprouting presence, height of sprouts
# Is there shrub/understory height data?
#   List of plots is attached as a .csv

lp <- read.csv(file.choose())
lp
# which lp plots were surveyed in 2020
lp$Plot[which(lp$Plot %in% S20plots)]

head(pAll)
dim(pAll)
pLS <- pAll[which(pAll$Plot.Orig %in% lp$Plot),]
dim(pLS)
unique(pLS$Species)
names(pLS)
names(pLS)[grep('fate',names(pLS))]

pLSx <- pLS[,c('Plot','Species','Num','Quad','X','Y','FirstSurvey','LastSurvey','LastType','Live.17','fate4.18','fate4.19','fate4.20','DBH_cm.17','DBH_cm.18','DBH_cm.19','DBH_cm.20','SA.BD_cm.17','SA.BD_cm.18','SA.BD_cm.19','SA.BD_cm.20','Basal.Resprout.Height_cm.18','Basal.Resprout.Height_cm.19','Basal.Resprout.Height_cm.20','Basal.Resprout.Count.18','Basal.Resprout.Count.19','Basal.Resprout.Count.20')]
pLSx <- pLSx[order(pLSx$Plot,pLSx$Quad,pLSx$Num),]
write.csv(pLSx,paste(od,'ForLauren.csv',sep='/'))

#########

# op.reset <- par(mfcol=c(1,1),mar=c(5,5,3,3))
# #fsCols <- c('blue','brown','orange','red')
# fsCols <- c('#00FF01','#FFB000','#FE6100','#DC267F')
# fsCols3 <- c("0.U" = '#00FF01',
#              "12.LM" = '#FE8800',
#              "3.H" = '#DC267F')
# vtShapes <- data.frame(vt=c('CON','MH','Mix3','WO'),shp=c(18,15,17,16))

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

# Read in species codes - in this script, called 'Species'
spAtt <- read.csv("https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/SpeciesData/all-spp-attributes.csv")

dim(spAtt)
head(spAtt)
tail(spAtt)
table(spAtt$FuncGroup)



######
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

table(tAll$Plot)
table(tAll$Survey)
ptAll <- tAll[which(tAll$Plot=='PPW1340' & tAll$Survey=='Plot'),]
dim(ptAll)
names(ptAll)[grep('Tag',names(ptAll))]
write.csv(file='ptAll[1:6,c('Plot','Quad.13','Type.13','X_cm.13','Y_cm.13','Num','Species','Live.17','fate3.18','fate3.19','fate3.20','BABH.17','BAd10.17','Basal.Resprout.Count.18','Basal.Resprout.Height_cm.18','Tag.Pulled.18','Tag.Pulled.19','Tag.Pulled.20')]


