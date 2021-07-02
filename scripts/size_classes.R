# make demographic distributions for species
rm(list=ls())
source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')

source('scripts/PW_functions_local.R')

seju.data<-get.seju.data(2018)
sapling<-plants.by.plot(year=2018,type="SA")
indv.data <- get.indv.data(year = 2018)


tree.data <- indv.data[which(indv.data$Type=='TR'),]

head(seju.data)
head(sapling)
head(tree.data)
tree.data$dbh <- 2*sqrt(tree.data$Basal.Area/pi)
max(tree.data$dbh,na.rm=T)

count.by.size <- function(id=tree.data,species=select.species,sr=c(1,10)) {
  ssel <- which(id$Species %in% species)
  dsel <- which(id$dbh>=sr[1] & id$dbh<sr[2])
  rsel <- intersect(ssel,dsel)
  return(length(rsel))
}

cuts <- c(1,2,4,8,16,32,64,128)

all.spp <- sort(unique(tree.data$Species))
all.spp

# if one species, enter code in first line and figure legend name in second line. select.species can also take lists of species, including all.spp to get all data
select.species <- 'QUEDOU'
sp.name <- 'BlueOak'
{
  sc <- data.frame(class=c('seedling','juvenile','sapling','T1','T2','T3','T4','T5','T6','T7'),species=sp.name,count=NA)
  sc$count[1] <- sum(seju.data$Num.Seedlings[which(seju.data$Species %in% select.species)])
  sc$count[2] <- sum(seju.data$Num.Juveniles[which(seju.data$Species %in% select.species)])
  sc$count[3] <- sum(sapling$Count[which(sapling$Species %in% select.species)])
  sc
  
  nclass <- length(cuts)-1
  i=1
  for (i in 1:nclass) {
    bks <- cuts[i:(i+1)]
    xx <- count.by.size(sr=bks)
    sc$count[i+3] <- xx
  }
  sc
}
png(paste('figures/sizedist_2018_',sp.name,'.png',sep=''),width = 800,height = 600)
if (select.species=='UMBCAL') barplot(sc$count,main=sp.name,ylim=c(0,1000)) else barplot(sc$count,main=sp.name)
dev.off()


## compare years
id18 <- indv.data
id13 <- get.indv.data(year = 2013)

id13[id13$Num==2976,]
id18[id18$Num==2976,]
