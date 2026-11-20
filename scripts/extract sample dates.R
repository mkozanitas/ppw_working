rm(list=ls())
require(RCurl)
require('lubridate')
get.plot<-function(year=2013){ 
  # list of file names
  if(year>=2018) plot.list<- paste("PPW",c(1301:1350,1851:1854), sep="")
  else plot.list<- paste("PPW",1301:1350, sep="")
  return(plot.list)
}

prefix <- 'https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/'
options(stringsAsFactors=FALSE) 

#getURL('https://github.com/dackerly/PepperwoodVegPlots/tree/master/2020/Woody2020/Data/OriginalCSV/Woody')
#dir.list <- dir(paste(prefix,year,"/Woody",year,"/Data/OriginalCSV/Woody/",sep=''))
year=2013
file.list <-(paste(prefix,year,"/Woody",year,"/Data/OriginalCSV/Woody/WoodySurvey",year,"_", sep='')) 

plot.list<-get.plot(year=year)

infile <- paste(file.list[1], plot.list[1], ".csv", sep='')
tmp <- read.csv(text=getURL(infile, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), na.strings=c("","NA"),nrows=1,header=F)
tmp
tmp[1,][grep('2013',tmp[1,])]

mega.data<-lapply(paste(file.list, plot.list, ".csv", sep=''), function(x) read.csv(text=getURL(x, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), na.strings=c("","NA") , skip=3)) 
names(mega.data) <- plot.list # assigns a name to each item of the list

first.line<-lapply(paste(file.list, plot.list, ".csv", sep=''), function(x) read.csv(text=getURL(x, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), na.strings=c("","NA") , nrow=1,header=F)) 
names(first.line) <- plot.list # assigns a name to each item of the list

sample.date <- NA
for (i in 1:length(first.line)) sample.date[i] <- first.line[[i]][11]
sample.date <- do.call(rbind,sample.date)

mYear <- substr(sample.date,1,4)
mDate <- round(decimal_date(mdy(paste(substr(sample.date,5,6),
                       substr(sample.date,7,8)
                       ,mYear,sep='/'))),2)
cbind(mYear,mDate)



