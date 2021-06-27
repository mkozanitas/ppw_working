# make demographic distributions for species
rm(list=ls())
source('scripts/PWFunctions_load.R')
#source('scripts/PWfunctions_GitHub_local.R')

source('scripts/PW_functions_local.R')

indv.data <- get.indv.data(year = 2013,branches=T)
head(indv.data)

sap.data <- indv.data[indv.data$Type=='SA',]
head(sap.data)
tail(sap.data)
dim(sap.data)

sapNum <- strsplit(as.character(sap.data$Num),split = '.',fixed=T)

sap.data$mainStem <- TRUE
head(sap.data)
for (i in 1:length(sapNum)) if (length(sapNum[[i]]) > 1) sap.data$mainStem[i] <- F


plot(SA.Stump.Height_cm~SA.Stump.BD_cm,data=sap.data[sap.data$mainStem==T,],xlim=c(0,5))
fit <- lm(SA.Stump.Height_cm~SA.Stump.BD_cm,data=sap.data[sap.data$mainStem==T,])
abline(fit)

plot(SA.Stump.BD_cm~SA.Stump.Height_cm,data=sap.data[sap.data$mainStem==T,],ylim=c(0,5))
fit <- lm(SA.Stump.BD_cm~SA.Stump.Height_cm,data=sap.data[sap.data$mainStem==T,])
abline(fit)
summary(fit)
     
plot(log(SA.Stump.BD_cm)~log(SA.Stump.Height_cm),data=sap.data[sap.data$mainStem==T,])
fit <- lm(log(SA.Stump.BD_cm)~log(SA.Stump.Height_cm),data=sap.data[sap.data$mainStem==T,])
abline(fit)
summary(fit)
exp(5.5)

#species specific
ssel <- 'UMBCAL'
plot(SA.Stump.BD_cm~SA.Stump.Height_cm,data=sap.data[sap.data$mainStem==T & sap.data$Species==ssel,],ylim=c(0,5),main=ssel)
fit <- lm(SA.Stump.BD_cm~SA.Stump.Height_cm,data=sap.data[sap.data$mainStem==T & sap.data$Species==ssel,])
abline(fit)
summary(fit)
