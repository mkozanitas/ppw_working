# bark thickness
bt <- read.csv('input_data/barkthickness.csv')
dim(bt)
head(bt)
table(bt$species)

t10 <- tapply(bt$barkat10,bt$species,mean,na.rm=T)
t20 <- tapply(bt$barkat20,bt$species,mean,na.rm=T)
t30 <- tapply(bt$barkat30,bt$species,mean,na.rm=T)
t50 <- tapply(bt$barkat50,bt$species,mean,na.rm=T)

t10
names <- c('ARBMEN','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')
bt.out <- data.frame(names,t10,t20,t30,t50)
write.csv(bt.out,'input_data/barkthickness_spmeans.csv')
