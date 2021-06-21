## Script to extract seju data
rm(list=ls())
source('scripts/PWfunctions_GitHub_Tutorial.R')
source('scripts/PWfunctions_GitHub_local.R')

seedling.juvenile<-get.seju.data(year=2013)

#look at results
head(seedling.juvenile)
dim(seedling.juvenile)
sum(seedling.juvenile$Total.Number) 

write.csv(seedling.juvenile, "data/seju.data.csv", row.names=F)

#not sure how to break these up by quad, plants.by.plot() is an existing function 

tree<-plants.by.plot(year=2013, type="TR")
sapling<-plants.by.plot(year=2013,type="SA")
all<-plants.by.plot(year=2013,type="SA.TR")
babies<-plants.by.plot(year = 2013, type = "SEJU")


#from Megan's code...looks the same but totals remove quad distinction

seju.data<-get.seju.data(2013)
seju.data<-aggregate(cbind(Num.Seedlings, Num.Juveniles,Total.Number)~Species+Plot,   data=seju.data, FUN=sum)
seju.data<-seju.data[,c(2,1,3,4,5)]







