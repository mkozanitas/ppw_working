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
print('L19d refers to plants that were last seen in 2019 and the plot WAS surveyed in 2020')
print('L19l refers to plants that were last seen in 2019 and the plot WAS NOT surveyed in 2020')
pAll$LastSurvey <- NA
pAll$LastSurvey[which(!is.na(pAll$fate3.20))] <- 'L20'
pAll$LastSurvey[which(is.na(pAll$LastSurvey) & !is.na(pAll$fate3.19) & (pAll$Plot %in% S20plots))] <- 'L19d'
pAll$LastSurvey[which(is.na(pAll$LastSurvey) & !is.na(pAll$fate3.19) & (!pAll$Plot %in% S20plots))] <- 'L19l'
pAll$LastSurvey[which(is.na(pAll$LastSurvey) & !is.na(pAll$fate3.18))] <- 'L18'

print("table(pAll$LastSurvey,useNA='a')")
table(pAll$LastSurvey,useNA='a')

print("table(pAll$FirstSurvey,pAll$LastSurvey,useNA='a')")
table(pAll$FirstSurvey,pAll$LastSurvey,useNA='a')

#table(pAll$Species,pAll$LastSurvey)

pAll$SampGroup <- 'No.SAMP'
pAll$SampGroup[which(pAll$FuncGroup %in% c('EHRO','NS.Con','WHTO'))] <- 'Yes.SAMP'
pAll$SampGroup[which(pAll$Species=='HETARB')] <- 'Yes.SAMP'
pAll$SampGroup[which(pAll$Species=='TORCAL')] <- 'No.SAMP'

print("table(pAll$Func,pAll$LastSurvey,pAll$SampGroup)")
table(pAll$Func,pAll$LastSurvey,pAll$SampGroup)

print("Number to sample by plot")
table(pAll$Plot[which(pAll$SampGroup=='Yes.SAMP')],pAll$FuncGroup[which(pAll$SampGroup=='Yes.SAMP')])

# How many total plants to sample
print("total number to sample in 2026")
length(which(pAll$LastSurvey %in% c('L19l','L20') & (pAll$SampGroup=='Yes.SAMP')))

sink()

# trial data sheets for 3 plots
# DBH (tree) or d10 (sap) when last seen
# fate when last seen

names(pAll)
pAll$LastType <- NA
pAll$LastDiam <- NA
pAll$LastFate3 <- NA
pAll$X <- NA
pAll$Y <- NA

r20 <- which(pAll$LastSurvey=='L20')
r19 <- which(pAll$LastSurvey %in% c('L19d','L19l'))
r18 <- which(pAll$LastSurvey=='L18')
pAll[r20,c('LastType','LastDiam','LastFate3')] <- pAll[r20,c('Type.20','DBH_cm.20','fate3.20')]
pAll[r19,c('LastType','LastDiam','LastFate3')] <- pAll[r19,c('Type.19','DBH_cm.19','fate3.19')]
pAll[r18,c('LastType','LastDiam','LastFate3')] <- pAll[r18,c('Type.18','DBH_cm.18','fate3.18')]
pAll$LastLive <- 'L'
pAll$LastLive[which(pAll$LastFate=='DN')] <- 'X'

r20 <- which(pAll$FirstSurvey=='F20')
r19 <- which(pAll$FirstSurvey=='F19')
r18 <- which(pAll$FirstSurvey=='F18')
r13 <- which(pAll$FirstSurvey=='F13')
pAll[r13,c('Quad','X','Y')] <- pAll[r13,c('Quad.13','X_cm.13','Y_cm.13')]
pAll[r18,c('Quad','X','Y')] <- pAll[r18,c('Quad.18','X_cm.18','Y_cm.18')]
pAll[r19,c('Quad','X','Y')] <- pAll[r19,c('Quad.19','X_cm.19','Y_cm.19')]
pAll[r20,c('Quad','X','Y')] <- pAll[r20,c('Quad.20','X_cm.20','Y_cm.20')]

pAll$TP <- 0
pAll$TP[which(pAll$Tag.Pulled.18!=0 | pAll$Tag.Pulled.19!=0 | pAll$Tag.Pulled.19!=0)] <- 1
table(pAll$TP)

sRows <- which(pAll$SampGroup=='Yes.SAMP' & pAll$TP==0)
length(sRows)
Sdat <- pAll[sRows,c('Plot','Num','Species','Quad','X','Y','FirstSurvey','LastSurvey','LastType','LastDiam','LastFate3','LastLive')]
head(Sdat)

pRows <- which(Sdat$Plot=='PPW1340')
write.csv(Sdat[pRows,],paste(od,'PPW1340.csv',sep='/'))

pRows <- which(Sdat$Plot=='PPW1335')
write.csv(Sdat[pRows,],paste(od,'PPW1335.csv',sep='/'))

pRows <- which(Sdat$Plot=='PPW1345')
write.csv(Sdat[pRows,],paste(od,'PPW1345.csv',sep='/'))

# 



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
