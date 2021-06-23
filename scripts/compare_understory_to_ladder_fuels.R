# compare seju counts to ladder fuels
rm(list=ls())


#read in pl data frame 

pl <- read.csv("data/pl_data.csv")
pl$DF.rem <- 0
pl$DF.rem[c(20,21,49,40,41,6,25)] <- 1 #from a random word document*

plot(cl~sapling.BA,data=pl[-28,],main = "Saplings (>50cm, <1cm DBH) by Basal Area",ylab='Classified ladder fuels',xlab='Sapling basal area (cm2)')
fit <- lm(cl~sapling.BA,data=pl[-28,])
abline(fit)
summary(fit)

plot(rl~sapling.BA,data=pl[-28,],main = "Saplings (>50cm, <1cm DBH) by Basal Area",pch=19)
points(rl~sapling.BA,data=pl[pl$DF.rem==1,],pch=19,col='red',cex=2)
fit <- lm(rl~sapling.BA,data=pl[-28,])
abline(fit)
summary(fit)

plot(rl~consap.BA,data=pl[-28,],main = "Saplings (>50cm, <1cm DBH) by Basal Area",pch=19)
points(rl~consap.BA,data=pl[pl$DF.rem==1,],pch=19,col='red',cex=2)
fit <- lm(rl~sapling.BA,data=pl[-28,])
abline(fit)
summary(fit)

plot(cl~sapling.count,data=pl[-28,],main = "Saplings (>50cm, <1cm DBH) by count totals")
fit <- lm(cl~sapling.count,data=pl[-28,])
abline(fit)
summary(fit)

plot(rl~sapling.count,data=pl[-28,],main = "Saplings (>50cm, <1cm DBH) by count totals")
fit <- lm(rl~sapling.count,data=pl[-28,])
abline(fit)
summary(fit)

plot(cl~TBA_1_3,data=pl, main = "Small trees (1-3cm DBH) by Basal Area")
fit <- lm(cl~TBA_1_3,data=pl)
abline(fit)
summary(fit)

plot(cl~TBA_3_5,data=pl, main = "Small trees (3-5cm DBH) by Basal Area")
fit <- lm(cl~TBA_3_5,data=pl)
abline(fit)
summary(fit)

plot(cl~TBA_5_7,data=pl, main = "Small trees (5-7cm DBH) by Basal Area")
fit <- lm(cl~TBA_5_7,data=pl)
abline(fit)
summary(fit)

plot(cl~TBA_7_10,data=pl, main = "Small trees (7-10cm DBH) by Basal Area")
fit <- lm(cl~TBA_7_10,data=pl)
abline(fit)
summary(fit)

plot(rl~TBA_1_3,data=pl, main = "Small trees (1-3cm DBH) by Basal Area")
fit <- lm(rl~TBA_1_3,data=pl)
abline(fit)
summary(fit)

plot(rl~TBA_3_5,data=pl, main = "Small trees (3-5cm DBH) by Basal Area")
fit <- lm(rl~TBA_3_5,data=pl)
abline(fit)
summary(fit)

plot(rl~TBA_5_7,data=pl, main = "Small trees (5-7cm DBH) by Basal Area")
fit <- lm(rl~TBA_5_7,data=pl)
abline(fit)
summary(fit)

plot(rl~TBA_7_10,data=pl, main = "Small trees (7-10cm DBH) by Basal Area")
fit <- lm(rl~TBA_7_10,data=pl)
abline(fit)
summary(fit)


# Multiple regression of contribution of saplings vs. small trees
head(pl)
cor(pl[,c('sapling.BA','TBA_1_3','TBA_3_5','TBA_5_7','TBA_7_10')])
pl$TBA_1_5 <- pl$TBA_1_3 + pl$TBA_3_5
pl$TBA_5_10 <- pl$TBA_5_7 + pl$TBA_7_10

cor(pl[-28,c('sapling.BA','TBA_1_5','TBA_5_10')])

fit <- glm(rl~sapling.BA + TBA_1_5,data=pl[-28,])
summary(fit)

plot(rl~sapling.BA,data=pl[-28,])
fit <- lm(rl~sapling.BA,data=pl[-28,])
summary(fit)

plot(rl~TBA_1_5,data=pl[-28,])
fit <- lm(rl~TBA_1_5,data=pl[-28,])
summary(fit)

plot(rl~TBA_1_3,data=pl[-28,])
fit <- lm(rl~TBA_1_3,data=pl[-28,])
summary(fit)

plot(rl~TBA_3_5,data=pl[-28,])
fit <- lm(rl~TBA_3_5,data=pl[-28,])
summary(fit)

#####This was used before creating main "pl" dataframe for pepperwood ladderfuels####

# read in seju data summarized by plot
seju.data.p <- read.csv('data/seju_data_p.csv')
dim(seju.data.p)

# read in ladder fuel data
ld <- read.csv('data/mean_ladder_fuel_54.csv')
dim(ld)

seju.data.p$ld <- ld$ld[1:50]
seju.data.p$rl <- ld$rl[1:50]
head(seju.data.p)

#plot(seju.data.p$ld~seju.data.p$Total.Number)
plot(ld~Total.Number,data=seju.data.p)
fit <- lm(ld~Total.Number,data=seju.data.p)
abline(fit)

# compare saplings
saplings.p <- read.csv('data/saplings_p.csv')
dim(saplings.p)

saplings.p$ld <- ld$ld[1:50]
saplings.p$rl <- ld$rl[1:50]
head(saplings.p)

# run these lines to replace NAs with zeros
saplings.p$Basal.Area[which(is.na(saplings.p$Basal.Area))] <- 0
saplings.p$Count[which(is.na(saplings.p$Count))] <- 0

plot(ld~Basal.Area,data=saplings.p[-28,],main = "Saplings (>50cm, <1cm DBH) by Basal Area")
fit <- lm(ld~Basal.Area,data=saplings.p[-28,])
abline(fit)
summary(fit)

plot(rl~Basal.Area,data=saplings.p[-28,],main = "Saplings (>50cm, <1cm DBH) by Basal Area")
fit <- lm(rl~Basal.Area,data=saplings.p[-28,])
abline(fit)
summary(fit)

plot(ld~Count,data=saplings.p[-28,],main = "Saplings (>50cm, <1cm DBH) by count totals")
fit <- lm(ld~Count,data=saplings.p[-28,])
abline(fit)
summary(fit)

plot(rl~Count,data=saplings.p[-28,],main = "Saplings (>50cm, <1cm DBH) by count totals")
fit <- lm(rl~Count,data=saplings.p[-28,])
abline(fit)
summary(fit)

# compare small trees
smalltrees.p <- read.csv('data/trees_1_5_p.csv')
dim(smalltrees.p)

smalltrees.p$ld <- ld$ld[1:50]
smalltrees.p$rl <- ld$rl[1:50]
head(smalltrees.p)

# run these lines to replace NAs with zeros
smalltrees.p$Basal.Area[which(is.na(smalltrees.p$Basal.Area))] <- 0
#smalltrees.p$Count[which(is.na(smalltrees.p$Count))] <- 0

smalltrees.p$log10.Basal.Area <- log10(smalltrees.p$Basal.Area+1)
smalltrees.p$log10.rl <- log10(smalltrees.p$rl+1)

plot(ld~Basal.Area,data=smalltrees.p, main = "Small trees (1-5cm DBH) by Basal Area")
fit <- lm(ld~Basal.Area,data=smalltrees.p)
abline(fit)
summary(fit)

#uncomment to identify specific points
#identify(smalltrees.p$ld~smalltrees.p$Basal.Area)

plot(rl~Basal.Area,data=smalltrees.p, main = "Small trees (1-5cm DBH) by Basal Area")
fit <- lm(rl~Basal.Area,data=smalltrees.p)
abline(fit)
summary(fit)

plot(log10.rl~log10.Basal.Area,data=smalltrees.p,main = "Small trees (1-5cm DBH) by Basal Area")
fit <- lm(log10.rl~log10.Basal.Area,data=smalltrees.p)
abline(fit)
summary(fit)



# raw vs. classified
plot(ld~rl,data=smalltrees.p)

## add together small trees and saplings
head(saplings.p)
head(smalltrees.p)
smalltrees.p$sap.smtr.BA <- saplings.p$Basal.Area+smalltrees.p$Basal.Area

plot(ld~sap.smtr.BA,data=smalltrees.p[-28,])
fit <- lm(ld~sap.smtr.BA,data=smalltrees.p[-28,])
abline(fit)
summary(fit)
identify(smalltrees.p$ld[-28]~smalltrees.p$sap.smtr.BA[-28])
