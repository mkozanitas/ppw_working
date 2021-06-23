## Script to extract plot data for plants of different sizes
rm(list=ls())
source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')

source('scripts/PW_functions_local.R')

#Extract and summarize saplings
sapling<-plants.by.plot(year=2013,type="SA")
#all<-plants.by.plot(year=2013,type="SA.TR")

# make summary files of saplings
head(sapling)
sapling.p <- aggregate(cbind(Basal.Area, Count)~Plot, data=sapling, FUN=sum)
head(sapling.p)
dim(sapling.p)

plot.list <- get.plot(2013)

match.plots <- match(plot.list,sapling.p$Plot)
sapling.p <- cbind(plot.list,sapling.p[match.plots,])
sapling.p <- sapling.p[,-2]
names(sapling.p)[1] <- 'Plot'
rownames(sapling.p) <- 1:50
head(sapling.p)
write.csv(sapling.p,'data/saplings_p.csv')

# Extract and summarize trees
tree<-plants.by.plot(year=2013, type="TR")
head(tree)

## extract small trees
indv.data <- get.indv.data(year = 2013)
head(indv.data)
dim(indv.data)

diam.range <- c(1,5) #enter min and max diam of interest
ba.range <- pi*(diam.range/2)^2

small.trees <- indv.data[which(indv.data$Basal.Area>=ba.range[1] & indv.data$Basal.Area<ba.range[2]),]
dim(small.trees)

small.trees.p <- aggregate(Basal.Area~Plot,data=small.trees,FUN=sum)
head(small.trees.p)

plot.list <- get.plot(2013)

match.plots <- match(plot.list,small.trees.p$Plot)
small.trees.p <- cbind(plot.list,small.trees.p[match.plots,])
small.trees.p <- small.trees.p[,-2]
names(small.trees.p)[1] <- 'Plot'
rownames(small.trees.p) <- 1:50
head(small.trees.p)

out.file <- paste('data/trees_',diam.range[1],'_',diam.range[2],'_p.csv',sep='')
write.csv(small.trees.p,out.file)

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


# %in% function
which(indv.data$Species %in% c('PSEMEN'))

#notes for buikding a data frame 
x=data.frame(c1=1:10)
x
x$c2 <- 11:20
x
x$c3 <- NA
x
x$c4 <- c(1,2,3)
x

ladderfuels=data.frame(c1=1:50)
ladderfuels$c2 <- sapling
