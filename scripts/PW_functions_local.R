### local function versions for trouble shooting
plants.by.plot <- function(year, type){
  plot.list<- get.plot()
  if(type=="SEJU"){
    seju<-get.seju.data(year=year)
    seju<-with(seju, aggregate(cbind(Num.Seedlings,Num.Juveniles,Total.Number)~Species+Plot+Year, FUN=sum))        
    seju<-seju[,c(2,1,3:6)]
    seju$Year<-year
    return(seju)                 
  }
  else{
    indv.data<-suppressWarnings(get.indv.data(year=year))
    if(type!="SA.TR") {
      count<-with(indv.data[indv.data$Type==type,], aggregate(Num~Species+Plot+Type, FUN=function(x)  length(unique(x))))
      area<-with(indv.data[indv.data$Type==type,], aggregate(Basal.Area~Species+Plot+Type, FUN=sum))
      output<-merge(area,count, by=c("Plot","Species","Type"))
      output$Year<-year
      colnames(output)<-c("Plot", "Species", "Type", "Basal.Area","Count","Year")}
    else{
      count.sap<-with(indv.data[indv.data$Type=="SA",], aggregate(Num~Species+Plot+Type, FUN=function(x)  length(unique(x))))
      area.sap<-with(indv.data[indv.data$Type=="SA",], aggregate(Basal.Area~Species+Plot+Type, FUN=sum))
      output.sap<-merge(area.sap,count.sap,by=c("Plot","Species","Type"))
      count.tr<-with(indv.data[indv.data$Type=="TR",], aggregate(Num~Species+Plot+Type, FUN=function(x)  length(unique(x))))
      area.tr<-with(indv.data[indv.data$Type=="TR",], aggregate(Basal.Area~Species+Plot+Type, FUN=sum))
      output.tr<-merge(area.tr,count.tr,by=c("Plot","Species","Type"))
      output<-merge(output.sap,output.tr,by=c("Plot","Species"),all=T)
      output<-output[-c(3,6)]
      output$Year<-year
      colnames(output)<-c("Plot","Species","SA.Basal.Area","SA.Count","TR.Basal.Area","TR.Count","Year")
      output[is.na(output)] <- 0
    }
    return(output)
  }   
}

get.indv.data <- function(year, stump=F, orig.dead=F, survival=F, bsprout=F, epicormic=F, apical=F, canopy=F, bsprout.height=F, bsprout.count=F, tag.pulled=F, keep.999=F, branches=F, prefix='https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/',plot.list=NA){
  options(stringsAsFactors=FALSE)  
  getURL('https://github.com/dackerly/PepperwoodVegPlots/tree/master/2020/Woody2020/Data/OriginalCSV/Woody')
  dir.list <- dir(paste(prefix,year,"/Woody",year,"/Data/OriginalCSV/Woody/",sep=''))
  file.list <-(paste(prefix,year,"/Woody",year,"/Data/OriginalCSV/Woody/WoodySurvey",year,"_", sep='')) 
  if (is.na(plot.list)) plot.list<-get.plot(year=year)
  
  mega.data<-lapply(paste(file.list, plot.list, ".csv", sep=''), function(x) read.csv(text=getURL(x, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), na.strings=c("","NA") , skip=3)) 
  names(mega.data) <- plot.list 
  
  i=1
  for (i in 1:length(mega.data)) 
  {
    Plot<-plot.list[i]
    mega.data[[i]]<-cbind(Plot=Plot, mega.data[[i]]) 
    
    if(year>=2018) {
      mega.data[[i]] <- mega.data[[i]][,c(1:4,7,13:16,5:6,20,8:12,17:19)]
      colnames(mega.data[[i]])<-c(
        "Plot", "Quad", "Type", "Num", "Species",
        "SA.Stump.Height_cm","SA.Stump.BD_cm",
        "SA.Branch.Num","DBH_cm", "X_cm", "Y_cm",
        "Notes","Survival", "Basal.Resprout",
        "Epicormic.Resprout","Apical.Growth",
        "Canopy.Percent", "Basal.Resprout.Height_cm",
        "Basal.Resprout.Count", "Tag.Pulled")
      if(tag.pulled==F) mega.data[[i]] <- mega.data[[i]][,-20]
      if(bsprout.count==F) mega.data[[i]] <- mega.data[[i]][,-19]
      if(bsprout.height==F) mega.data[[i]] <- mega.data[[i]][,-18]
      if(canopy==F) mega.data[[i]] <- mega.data[[i]][,-17]
      if(apical==F) mega.data[[i]] <- mega.data[[i]][,-16]
      if(epicormic==F) mega.data[[i]] <- mega.data[[i]][,-15]
      if(bsprout==F) mega.data[[i]] <- mega.data[[i]][,-14]
      if(survival==F) mega.data[[i]] <- mega.data[[i]][,-13]    
      if(keep.999==F) mega.data[[i]] <- mega.data[[i]][mega.data[[i]]$Num>=1000,]
    } else {
      mega.data[[i]]<-mega.data[[i]][,c(1:5,7:14)] 
      colnames(mega.data[[i]])<-c("Plot", "Quad", "Type", "Num", "Species", "Dead.Stump", 
                                  "SA.Stump.Height_cm", "SA.Stump.BD_cm", "SA.Branch.Num", "DBH_cm", "X_cm", "Y_cm", "Notes") 
      # subset stumps and original dead
      if(stump==F & orig.dead==T) mega.data[[i]]<-subset(mega.data[[i]], subset=(mega.data[[i]]$Dead.Stump=="D" | is.na(mega.data[[i]]$Dead.Stump))) 
      # subset original dead individuals
      if(stump==T & orig.dead==F) mega.data[[i]]<-subset(mega.data[[i]], subset=(mega.data[[i]]$Dead.Stump=="S" | is.na(mega.data[[i]]$Dead.Stump)))
      #subset both 
      if(stump==F & orig.dead==F) {mega.data[[i]]<-subset(mega.data[[i]], subset=(is.na(mega.data[[i]]$Dead.Stump)))
      #mega.data[[i]]<-mega.data[[i]][,-6]}
      }
    }
  }
  
  # turn list of dataframes into a single large dataframe called indv.data for easier manipulations
  indv.data<-do.call(rbind, mega.data)
  # get rid of extra rows
  indv.data<-indv.data[!is.na(indv.data$Type),]
  # change SA.Stump.BD_cm into numeric
  indv.data$SA.Stump.BD_cm <- suppressWarnings(as.numeric(indv.data$SA.Stump.BD_cm))
  # make indv.data$TreeNum numeric
  indv.data$Num<-as.numeric(indv.data$Num)
  # cleaning up Types 
  indv.data[indv.data$Type==" TR", "Type"]<-"TR"
  indv.data[indv.data$Type=="SA ","Type"]<-"SA"
  indv.data[indv.data$Type=="SA  ","Type"]<-"SA"
  indv.data[indv.data$Type=="S","Type"]<-"SA"
  indv.data[indv.data$Type=="AS","Type"]<-"SA"
  
  # fix up Species names for bad data entry, etc.
  indv.data[which(indv.data$Species=="CEOCUN"), "Species"]<-"CEACUN"
  indv.data[which(indv.data$Species=="UNKN27"), "Species"]<-"CEACUN"
  indv.data[which(indv.data$Species=="ARBMEN "), "Species"]<-"ARBMEN"
  indv.data[which(indv.data$Species=="ARCMAN "), "Species"]<-"ARCMAN"
  indv.data[which(indv.data$Species=="CERBET "), "Species"]<-"CERBET"
  indv.data[which(indv.data$Species=="QUEAGR "), "Species"]<-"QUEAGR"
  indv.data[which(indv.data$Species=="QUEAGRI"), "Species"]<-"QUEAGR"
  indv.data[which(indv.data$Species=="QUEARG"), "Species"]<-"QUEAGR"
  indv.data[which(indv.data$Species=="QUEGAR "), "Species"]<-"QUEGAR"
  indv.data[which(indv.data$Species=="QUEKEL "), "Species"]<-"QUEKEL"
  indv.data[which(indv.data$Species=="PIEMEN"), "Species"]<-"PSEMEN"
  indv.data[which(indv.data$Species=="PSMEN"), "Species"]<-"PSEMEN"
  indv.data[which(indv.data$Species=="UMBCAL "), "Species"]<-"UMBCAL"
  indv.data[which(indv.data$Species=="UMCAL"), "Species"]<-"UMBCAL"
  indv.data[which(indv.data$Species=="UNKN30"), "Species"]<-"UNK"
  indv.data[which(indv.data$Species=="QUEDEC"), "Species"]<-"QUEdec"
  
  # Change individuals originally identified as "QUEDEC" to species-level indentification 
  if(year==2013) {
    AUG.ID<-read.csv(text=getURL(paste(prefix, "2013/OakID2013/AUG_Species.csv", sep=''), followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))
    for(i in 1:dim(AUG.ID)[1]){
      indv.data[indv.data$Num %in% AUG.ID$Num[i], "Species"] <-AUG.ID$Species[i]
    }
  }
  
  # change coordinates to all be positive
  indv.data$Y_cm<-suppressWarnings(as.numeric(indv.data$Y_cm))
  indv.data[which(indv.data$X_cm<0),"X_cm"]<-indv.data[which(indv.data$X_cm<0),"X_cm"]+500
  indv.data[which(indv.data$Y_cm<0),"Y_cm"]<-indv.data[which(indv.data$Y_cm<0),"Y_cm"]+500
  
  # condense each individual into a single row 
  if (branches==F){
    indv.data$Num<-floor(indv.data$Num)
    if(year==2018) {
      indv.data[indv.data$Type=="SA" & is.na(indv.data$SA.Stem.Num), "SA.Stem.Num" ] <- 0
      indv.data[indv.data$Type=="SA", "SA.Stem.Num" ]<- indv.data[indv.data$Type=="SA", "SA.Stem.Num" ] + 1
    }
    indv.data[indv.data$Type=="SA", "Basal.Area" ]<-(((indv.data[indv.data$Type=="SA", "SA.Stump.BD_cm"]/2)^2)*(pi))*(indv.data[indv.data$Type=="SA","SA.Branch.Num"])
    
    indv.data[indv.data$Type=="TR", "Basal.Area"]<- (((indv.data[indv.data$Type=="TR","DBH_cm"]/2)^2)*(pi))
    # get rid of TS rows
    indv.data<-indv.data[indv.data$Type!="TS",]
    library(data.table)
    indv.data<-data.table(indv.data)
    indv.data[,Basal.Area:=sum(Basal.Area),by="Num"]
    if(year==2018) {
      if(survival==TRUE) indv.data[,Survival:=max(Survival),by="Num"]
      if(bsprout==TRUE) indv.data[,Basal.Resprout:=max(Basal.Resprout),by="Num"]
      if(epicormic==TRUE) indv.data[,Epicormic.Resprout:=max(Epicormic.Resprout),by="Num"]
      if(apical==TRUE) indv.data[,Apical.Growth:=max(Apical.Growth),by="Num"]
      if(canopy==TRUE) indv.data[,Canopy.Percent:=max(Canopy.Percent),by="Num"] #### BIT OF A FUDGE
    }
    # get rid of duplicated tag numbers (?)
    indv.data<-indv.data[!duplicated(indv.data$Num,incomparables=NA),] 
    
    # Get rid of the extra columns "DBH_cm", ... 
    indv.data<-as.data.frame(indv.data)
    indv.data<-indv.data[,c("Plot","Quad","Type","Num","Species","X_cm","Y_cm","Basal.Area")]
  } else {
    indv.data[indv.data$Type=="SA", "Basal.Area" ]<-(((indv.data[indv.data$Type=="SA", "SA.Stump.BD_cm"]/2)^2)*(pi))*(indv.data[indv.data$Type=="SA","SA.Branch.Num"])
    indv.data[indv.data$Type=="TR", "Basal.Area"]<- (((indv.data[indv.data$Type=="TR","DBH_cm"]/2)^2)*(pi))
  }
  
  if (year=="none"){
    indv.data<-subset(indv.data, subset=(indv.data$Dead.Stump=="D" | indv.data$Dead.Stum=="S"))  
  } else {
    if(year!=2012){
      indv.data$Year<-year
      #########
      # MORTALITY FUNCTIONS SHOULD NOT EVEN BE HERE SHOULD THEY?
      # MAKE A NEW FUNCTION get.indv.mortality.data()
      #    year<-2013:year
      #    dead<-kill.trees(year=year)
      #    indv.data<-indv.data[!(indv.data$Num %in% dead$Num),]
      #########
    }
  }
  
  row.names(indv.data) <- NULL
  
  # remove PSEMEN trees that were accidently removed by chainsaw crews in winter 2017 in plot PPW1307 
  indv.data<-indv.data[which(!(indv.data$Plot=="PPW1307" & indv.data$Num==1108 | indv.data$Num==1115 |indv.data$Num==1118 |indv.data$Num==1121 | indv.data$Num==1127 | indv.data$Num==1169 | indv.data$Num==1170 |indv.data$Num==1174)),]
  
  return(indv.data)
}
