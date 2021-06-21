## Script to extract seju data
rm(list=ls())
source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')




##### GIVING ERRORS BELOW....
#not sure how to break these up by quad, plants.by.plot() is an existing function 

source('scripts/PW_functions_local.R')
tree<-plants.by.plot(year=2013, type="TR")
sapling<-plants.by.plot(year=2013,type="SA")
all<-plants.by.plot(year=2013,type="SA.TR")

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





