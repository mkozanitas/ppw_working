## Script to extract seju data
rm(list=ls())
source('scripts/PWfunctions_GitHub_Tutorial.R')
source('scripts/PWfunctions_GitHub_local.R')

#from Megan's code...looks the same but totals remove quad distinction

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


##### GIVING ERRORS BELOW....
#not sure how to break these up by quad, plants.by.plot() is an existing function 

tree<-plants.by.plot(year=2013, type="TR")
sapling<-plants.by.plot(year=2013,type="SA")
all<-plants.by.plot(year=2013,type="SA.TR")
babies<-plants.by.plot(year = 2013, type = "SEJU")










