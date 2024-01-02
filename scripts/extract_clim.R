rm(list=ls())
library(raster)
## pull up plot info
plotclim <- read.csv('input_data/plot_info/clim_pts.csv')
head(plotclim)
tail(plotclim)

plotinfo <- read.csv('input_data/plot_info/plot.info.csv')
head(plotinfo)
tail(plotinfo)

names(plotclim)

cwd <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/bcm_outputs/aea/cwd1981_2010_ave.asc')
cwd[cwd<0] <- 0
plot(cwd,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$cwd <- extract(cwd,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/marrad.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$marrad <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/model3.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$model3 <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/PLP_100.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$plp100 <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/PLP_500.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$plp500 <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/topoidx.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$topoid <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_100.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$tpi100 <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_500.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$tpi500 <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_1000.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$tpi1000 <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/pwd_base/pw_10m_t2.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$dem <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/pwd_base/pwd_10m_t1.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$dem <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

prx <- raster('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/slpdegree.asc')
plot(prx,xlim=c(-241000,-230000),ylim=c(63000,72000))
points(plotinfo$AEA.x,plotinfo$AEA.y)
plotinfo$slope <- extract(prx,plotinfo[,c('AEA.x','AEA.y')])

write.csv(plotinfo,'input_data/plot_info/plot.info.csv')

names(plotinfo)
pairs(plotinfo[,9:19])
