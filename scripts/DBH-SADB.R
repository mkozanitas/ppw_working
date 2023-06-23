# file.choose gives you mac OPEN dialog to choose file
bd <- read.csv('data/BDtoDBH.csv')
dim(bd)
head(bd)
names(bd)
names(bd)[15] <- 'DBH.cm'
names(bd)[13] <- 'SA.BD.cm'

#create log10 vars for dbh and stem bd
bd$lDBH <- log10(bd$DBH.cm)
bd$lSABD <- log10(bd$SA.BD.cm)

# first look, log transformed, with regression line, 1:1 line, and threshold for DBH=1 (i.e. log DBH=0)
plot(lDBH~lSABD,data=bd,log='',xlim=c(-1,2),ylim=c(-1,2))
fit <- lm(lDBH~lSABD,data=bd)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(fit)

# first look, not logged
plot(DBH.cm~SA.BD.cm,data=bd,log='')
fit <- lm(DBH.cm~SA.BD.cm,data=bd)
abline(0,1,lty=2)
abline(h=1,lty=2)
abline(fit)

# remove points with DBH>SA.BD and the four points with really low DBH vs. SA.BD.cm
rm1 <- which(bd$DBH.cm>=bd$SA.BD.cm)
bd[rm1,]

rm2 <- which(bd$DBH.cm<10 & bd$SA.BD.cm>15)
bd[rm2,]

plot(lDBH~lSABD,data=bd)
points(lDBH~lSABD,data=bd[rm1,],pch=19)
points(lDBH~lSABD,data=bd[rm2,],pch=19)

# new data frame with those 13 points removed
bd2 <- bd[!1:nrow(bd) %in% c(rm1,rm2),]
dim(bd)
dim(bd2)

# plot log data again
plot(lDBH~lSABD,data=bd2,log='',xlim=c(-1,2),ylim=c(-1,2))
fit <- lm(lDBH~lSABD,data=bd)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(fit)

#linear data
plot(DBH.cm~SA.BD.cm,data=bd2,log='')
fit <- lm(DBH.cm~SA.BD.cm,data=bd)
abline(0,1,lty=2)
abline(fit)

#linear data focusing on points near 0
plot(DBH.cm~SA.BD.cm,data=bd2,log=''
,xlim=c(0,5),ylim=c(0,5))
fit <- lm(DBH.cm~SA.BD.cm,data=bd)
fit
abline(0,1,lty=2)
abline(h=1,lty=2)
abline(fit)

# log regression is better in terms of the distribution of points. But by definition logs can't go to zero. In contrast, for this data set, DBH  really does go to zero for saplings with height < 1.5m, so I think linear is better. So we would use the equation above (fit) to convert linear sapling basal area to dbh.
fit
# DBH = SA.BD * 0.64426 - 0.06439

#Looking at the slope, you can see the regression almost converges on 0, so we're not going to get dbh<0 (which would be hard to explain).

# next step is to read in all the sapling data, convert to dbh with the equation above. Then those converted values would be substituted for all of our size-based analyses

