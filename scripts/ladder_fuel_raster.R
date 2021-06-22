## reading in ladder fuel data
require('raster')
require('maptools')
require('sp')
require('rgeos')
require("rgdal")

## input 64m from Sonoma Veg and extract per 50 plots
ld <- raster('input_data/Classified_Ladder_Fuels/ladder.tif')
ld
plot(ld)
proj4string(ld)

PWD <- readShapeSpatial('input_data/PPshapefile-teale-albers/pepperwood')
plot(PWD)
proj4string(PWD) <- CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000')

PWDl <- spTransform(PWD,proj4string(ld))
plot(ld)
lines(PWDl)

Pld <- crop(ld,PWDl)
plot(Pld)
lines(PWDl)

p54 <- readShapeSpatial('input_data/shapefiles/vegplots-54-20m-aea/vegplots-54-20m')
proj4string(p54) <- proj4string(PWD)
p54l <- spTransform(p54,proj4string(ld))

plot(Pld)
plot(PWDl,bg='transparent',add=T)
plot(p54l,bg="transparent", add=TRUE)

ldf <- extract(Pld,p54l)
str(ldf)
ldf

ldf.sum <- data.frame(plot=1:54,ld=NA)
head(ldf.sum)
for (i in 1:54) ldf.sum$ld[i] <- mean(ldf[[i]])
ldf.sum



# input 20m provided by Ryan and extract for 20 m plots
## TURNS OUT TO BE THE SAME DATA - COMMENTED OUT
# ld2 <- raster('input_data/lf_clip/')
# 
# proj4string(ld2)
# proj4string(p54l)
# 
# plot(ld2)
# plot(p54l,bg="transparent", add=TRUE)
# 
# ldf2 <- extract(ld2,p54l)
# class(ldf2)
# ldf2[[1]]
# 
# ldf.sum$ld2 <- NA
# for (i in 1:54) ldf.sum$ld2[i] <- mean(ldf2[[i]])
# plot(ld2~ld,data=ldf.sum)


# RAW Ladder fuels data
rl <- raster('input_data/Raw_Ladder_Fuels/COUNTY_one_to_four.tif')
Prl <- crop(rl,PWDl)
plot(Prl)
plot(PWDl,bg='transparent',add=T)
plot(p54l,bg="transparent", add=TRUE)

rlf <- extract(Prl,p54l)
str(rlf)
rlf

ldf.sum$rl <- NA
head(ldf.sum)
for (i in 1:54) ldf.sum$rl[i] <- mean(rlf[[i]])
ldf.sum
plot(ldf.sum$ld,ldf.sum$rl)

write.csv(ldf.sum,'data/mean_ladder_fuel_54.csv')
