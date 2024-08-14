# file.choose gives you mac OPEN dialog to choose file
bd <- read.csv('input_data/BDtoDBH.csv')
dim(bd)
head(bd)
names(bd)
names(bd)[15] <- 'DBH.cm'
names(bd)[13] <- 'SA.BD.cm'
table(bd$Plot,bd$Type)

#create log10 vars for dbh and stem bd
bd$lDBH <- log10(bd$DBH.cm)
bd$lSABD <- log10(bd$SA.BD.cm)

# first look, log transformed, with regression line, 1:1 line, and threshold for DBH=1 (i.e. log DBH=0)
plot(lSABD~lDBH,data=bd,log='',xlim=c(-1,2),ylim=c(-1,2))
fit <- lm(lSABD~lDBH,data=bd)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(fit)

# first look, not logged
plot(SA.BD.cm~DBH.cm,data=bd,log='')
fit <- lm(SA.BD.cm~DBH.cm,data=bd)
abline(0,1,lty=2)
abline(h=1,lty=2)
abline(fit)

# remove points with DBH>SA.BD and the four points with really low DBH vs. SA.BD.cm
rm1 <- which(bd$DBH.cm>=bd$SA.BD.cm)
bd[rm1,]

rm2 <- which(bd$DBH.cm<10 & bd$SA.BD.cm>15)
bd[rm2,]

# saplings in this data set are a non-random sample where they seemed to be a tree, and had a dbh measured, and then realized they were too small. 
rm3 <- which(bd$Type=='SA')
bd[rm3,]
length(rm3)

plot(lSABD~lDBH,data=bd)
points(lSABD~lDBH,data=bd[rm1,],pch=19)
points(lSABD~lDBH,data=bd[rm2,],pch=19)
points(lSABD~lDBH,data=bd[rm3,],pch=19,col='red')

# new data frame with 4 outlier points removed. The 4 outliers are correctly entered from scan and taper from BD to DBH- dropping them
bd2 <- bd[!1:nrow(bd) %in% c(rm1),]
dim(bd)
dim(bd2)

# plot log data again
plot(lSABD~lDBH,data=bd2,log='',xlim=c(-1,2),ylim=c(-1,2))
fit <- lm(lSABD~lDBH,data=bd2)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(fit)
fit

#linear data
plot(SA.BD.cm~DBH.cm,data=bd2,log='')
fit.all <- lm(SA.BD.cm~DBH.cm,data=bd)
abline(0,1,lty=2)
abline(fit.all,col='red')

#linear data focusing on points near 0
plot(SA.BD.cm~DBH.cm,data=bd2,log='',xlim=c(0,5),ylim=c(0,5))
abline(fit.all,col='red')
fit.5 <- lm(SA.BD.cm~DBH.cm,data=bd2[which(bd2$SA.BD.cm<=5),])
abline(0,1,lty=2)
abline(h=1,lty=2)
abline(fit.5)

# log regression is better in terms of the distribution of points. But by definition logs can't go to zero. In contrast, for this data set, DBH really does go to zero for saplings with height < 1.5m, so I think linear is better. So we would use the equation above (fit) to convert linear sapling basal area to dbh.4
summary(fit.all)
length(complete.cases(bd[,c('SA.BD.cm','DBH.cm')]))
# D10 = DBH.cm * 1.176 + 1.070

# to revert back
# DBH.cm = (D10-1.070)/1.176

#Looking at the slope, you can see the regression almost converges on 0, so we're not going to get dbh<0 (which would be hard to explain).

# next step is to read in all the sapling data, convert to dbh with the equation above. Then those converted values would be substituted for all of our size-based analyses

