## Script to extract plot data for plants of different sizes
rm(list=ls())
source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')

source('scripts/PW_functions_local.R')

# Setup master dataframe for plot level analyses
plot.list <- get.plot(2013)
pl <- data.frame(order=1:50)
head(pl)
pl$order

pl$Plot=plot.list
head(pl)

ld <- read.csv('data/mean_ladder_fuel_54.csv')
dim(ld)
head(ld)

pl$cl <- ld$ld[1:50]
pl$rl <- ld$rl[1:50]
head(pl)

plot(pl$rl,pl$cl)

#Extract and summarize saplings
sapling<-plants.by.plot(year=2013,type="SA")
#all<-plants.by.plot(year=2013,type="SA.TR")

# make summary files of saplings
head(sapling)
sapling.p <- aggregate(cbind(Basal.Area, Count)~Plot, data=sapling, FUN=sum)
head(sapling.p)
dim(sapling.p)

match.plots <- match(pl$Plot,sapling.p$Plot)
sapling.p <- cbind(plot.list,sapling.p[match.plots,])
head(sapling.p)
sapling.p$Basal.Area[which(is.na(sapling.p$Basal.Area))] <- 0
sapling.p$Count[which(is.na(sapling.p$Count))] <- 0
head(sapling.p)

head(pl)
pl$sapling.count <- sapling.p$Count
pl$sapling.BA <- sapling.p$Basal.Area
head(pl)

# extract conifer counts and basal areas
head(sapling)
dim(sapling)
conifer.sapling <- sapling[which(sapling$Species %in% c('PSEMEN','TORCAL')),]
dim(conifer.sapling)

conifer.sapling.p <- aggregate(cbind(Basal.Area, Count)~Plot, data=conifer.sapling, FUN=sum)
head(conifer.sapling.p)
dim(conifer.sapling.p)

match.plots <- match(pl$Plot,conifer.sapling.p$Plot)
conifer.sapling.p <- cbind(plot.list,conifer.sapling.p[match.plots,])
head(conifer.sapling.p)
dim(conifer.sapling.p)

conifer.sapling.p$Basal.Area[which(is.na(conifer.sapling.p$Basal.Area))] <- 0
conifer.sapling.p$Count[which(is.na(conifer.sapling.p$Count))] <- 0
head(conifer.sapling.p)

head(pl)
pl$consap.count <- conifer.sapling.p$Count
pl$consap.BA <- conifer.sapling.p$Basal.Area
head(pl)

# extract shrub counts and basal areas
head(sapling)
dim(sapling)
table(sapling$Species)

# four unknown saplings to possible clean up later

# check shrub species list for future years
shrub.sapling <- sapling[which(sapling$Species %in% c('HETARB','CEACUN','FRACAL','ARCMAN','AMOCAL','BACPIL','RHACRO','SAMNIG','QUEBER','QUEBEGA','PRUILI','PRUCER','CORCOR','CEATHY','ADEFAS','CERBET','HOLDIS')),]
dim(shrub.sapling)
head(shrub.sapling)

shrub.sapling.p <- aggregate(cbind(Basal.Area, Count)~Plot, data=shrub.sapling, FUN=sum)
head(shrub.sapling.p)
dim(shrub.sapling.p)

match.plots <- match(pl$Plot,shrub.sapling.p$Plot)
shrub.sapling.p <- cbind(plot.list,shrub.sapling.p[match.plots,])
head(shrub.sapling.p)
dim(shrub.sapling.p)

shrub.sapling.p$Basal.Area[which(is.na(shrub.sapling.p$Basal.Area))] <- 0
shrub.sapling.p$Count[which(is.na(shrub.sapling.p$Count))] <- 0
head(shrub.sapling.p)

head(pl)
pl$shrubsap.count <- shrub.sapling.p$Count
pl$shrubsap.BA <- shrub.sapling.p$Basal.Area
head(pl)

pl$htsap.count <- pl$sapling.count - pl$consap.count - pl$shrubsap.count
pl$htsap.BA <- pl$sapling.BA - pl$consap.BA - pl$shrubsap.BA
head(pl)

## DONE WITH SAPLINGS!!

## MOVE ON TO TREES, by size class
indv.data <- get.indv.data(year = 2013)
head(indv.data)
dim(indv.data)

diam.range <- c(7,10) #enter min and max diam of interest
ba.range <- pi*(diam.range/2)^2

small.trees <- indv.data[which(indv.data$Basal.Area>=ba.range[1] & indv.data$Basal.Area<ba.range[2]),]
dim(small.trees)
table(small.trees$Species)

small.trees.p <- aggregate(Basal.Area~Plot,data=small.trees,FUN=sum)
head(small.trees.p)
dim(small.trees.p)

match.plots <- match(plot.list,small.trees.p$Plot)
small.trees.p <- cbind(plot.list,small.trees.p[match.plots,])
head(small.trees.p)

small.trees.p$Basal.Area[which(is.na(small.trees.p$Basal.Area))] <- 0
head(small.trees.p)

colname <- paste('TBA_',diam.range[1],'_',diam.range[2],sep='')
pl$x <- small.trees.p$Basal.Area
head(pl)
names(pl)[grep('x',names(pl))] <- colname
head(pl)

## pull out conifers
indv.data.archive <- indv.data


# subset to conifers
indv.data <- indv.data.archive[which(indv.data.archive$Species %in% c('PSEMEN','TORCAL')),]

diam.range <- c(7,10) #enter min and max diam of interest
ba.range <- pi*(diam.range/2)^2

small.trees <- indv.data[which(indv.data$Basal.Area>=ba.range[1] & indv.data$Basal.Area<ba.range[2]),]
dim(small.trees)
table(small.trees$Species)

small.trees.p <- aggregate(Basal.Area~Plot,data=small.trees,FUN=sum)
head(small.trees.p)
dim(small.trees.p)

match.plots <- match(plot.list,small.trees.p$Plot)
small.trees.p <- cbind(plot.list,small.trees.p[match.plots,])
head(small.trees.p)

small.trees.p$Basal.Area[which(is.na(small.trees.p$Basal.Area))] <- 0
head(small.trees.p)

colname <- paste('conTBA_',diam.range[1],'_',diam.range[2],sep='')
pl$x <- small.trees.p$Basal.Area
head(pl)
names(pl)[grep('x',names(pl))] <- colname
head(pl)

# subset to shrubs
indv.data <- indv.data.archive[which(indv.data.archive$Species %in% c('HETARB','CEACUN','FRACAL','ARCMAN','AMOCAL','BACPIL','RHACRO','SAMNIG','QUEBER','QUEBEGA','PRUILI','PRUCER','CORCOR','CEATHY','ADEFAS','CERBET','HOLDIS')),]
dim(indv.data)
head(indv.data)

diam.range <- c(1,3) #enter min and max diam of interest
ba.range <- pi*(diam.range/2)^2

small.trees <- indv.data[which(indv.data$Basal.Area>=ba.range[1] & indv.data$Basal.Area<ba.range[2]),]
dim(small.trees)
table(small.trees$Species)

small.trees.p <- aggregate(Basal.Area~Plot,data=small.trees,FUN=sum)
head(small.trees.p)
dim(small.trees.p)

match.plots <- match(plot.list,small.trees.p$Plot)
small.trees.p <- cbind(plot.list,small.trees.p[match.plots,])
head(small.trees.p)

small.trees.p$Basal.Area[which(is.na(small.trees.p$Basal.Area))] <- 0
head(small.trees.p)

colname <- paste('shrubTBA_',diam.range[1],'_',diam.range[2],sep='')
pl$x <- small.trees.p$Basal.Area
head(pl)
names(pl)[grep('x',names(pl))] <- colname
head(pl)
pl[,c(colname,'Plot')]

#### NOT USED FOR NOW
#Pull SEJU data and create summary files

seju.data<-get.seju.data(2013)

#look at results
head(seju.data)
dim(seju.data)
sum(seju.data$Total.Number) 

seju.data.sp<-aggregate(cbind(Num.Seedlings, Num.Juveniles,Total.Number)~Species+Plot, data=seju.data, FUN=sum)
head(seju.data.sp)
seju.data.sp<-seju.data.sp[,c(2,1,3,4,5)]
head(seju.data.sp)

seju.data.s<-aggregate(cbind(Num.Seedlings, Num.Juveniles,Total.Number)~Species, data=seju.data, FUN=sum)
seju.data.s

seju.data.p<-aggregate(cbind(Num.Seedlings, Num.Juveniles,Total.Number)~Plot, data=seju.data, FUN=sum)
head(seju.data.p)

write.csv(seju.data, "data/seju_data.csv", row.names=F)
write.csv(seju.data.sp, "data/seju_data_sp.csv", row.names=F)
write.csv(seju.data.p, "data/seju_data_p.csv", row.names=F)
write.csv(seju.data.s, "data/seju_data_s.csv", row.names=F)





