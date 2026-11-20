## This script takes the first steps in downloading and preparing the data used for the analysis published in Kozanitas et al. 2026 (DOI 10.1016/j.foreco.2026.123706). The data and some functions are stored on github. Note that during QA/QC of the data, sometimes we fix the data files, but often we make corrections in the script. Thus, it is important to run this script anytime an analysis is starting from the raw data. This script is fairly short and separate from the next steps because the queries to github are very slow.
# This script downloads two versions of the data - one at the individual tree level and the other at the stem level. Data from all years (2013-2020) is combined, and saved to new Robjects for future ues. The next steps in data analysis are conducted in K26_20_cleanData.R. 

##Note that the Kozanitas publication used the data at the stem level (i.e. branches=T)

rm(list=ls())

## LOAD PPW FUNCTIONS

# this code is from some delightful human on the internet who figured out how to source from GitHub; run this function in order to source in the script with all the Pepperwood Functions
library("RCurl")
source_https <- function(url, ...) { 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

# sources in all functions (described below) that allow access to the PPW Vegetation Plot Data

# once functions are uploaded to github, swap these two lines
source_https('https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/Analyses/PW_functions_GitHub_2026update.R')
########

## FOR NOW CODE TO EXPORT DATA WITH BRANCHES ALL COMMENTED OUT
# Extract plot data for 2013, 2018, 2019, 2020 and prepare for analysis

# Load the per-year data (without aggregating branches), then set to branches=F to collapse points
id13 <- get.indv.data(year = 2013,branches=F,stump=F,orig.dead=F)
dim(id13)
head(id13)
names(id13)
head(sort(id13$Num))

#to see all individual branches...points not collapsed use branches=T
id13b <- get.indv.data(year = 2013,branches=T,stump=F,orig.dead=F)
dim(id13b)
head(id13b)
head(sort(id13b$Num))
id13b[id13b$Num==3482,]

# sloppy coding to be able to turn a code snippet on and off
if (FALSE) 
{
  id13f <- get.indv.data(year = 2013, stump=T, orig.dead=T, survival=T, bsprout=T, epicormic=T, apical=T, canopy=T, bsprout.height=T, bsprout.count=T, tag.pulled=T, keep.999=T, branches=F)
  dim(id13f)
  head(id13f)
}

#indv.data.2013.B <- get.indv.data(year = 2013,branches=T)
#dim(indv.data.2013.B)
#head(indv.data.2013.B)

id18 <- get.indv.data(year = 2018,branches=F,keep.999 = T)
dim(id18)
head(id18)

id18b <- get.indv.data(year = 2018,branches=T,keep.999 = T)
dim(id18b)
head(id18b)

## run 50.troubleshoot interactively
####
id19 <- get.indv.data(year = 2019,branches=F)
dim(id19)
head(id19)

id19b <- get.indv.data(year = 2019,branches=T)
dim(id19b)
head(id19b)

plot.list.2020 <- get.plot(2020)
# now remove plots not sampled
plot.list.2020 <- plot.list.2020[-c(5,6,20,25,28,30,40:42,46,49,51:54)]
plot.list.2020

id20 <- get.indv.data(year = 2020,branches=F,plot.list=plot.list.2020)
dim(id20)
head(id20)

id20b <- get.indv.data(year = 2020,branches=T,plot.list=plot.list.2020)
dim(id20b)
head(id20b)

## write out all the data to an RData file, so it can be read and analyzed without going back to github
all.id <- list(id13,id18,id19,id20)
all.idb <- list(id13b,id18b,id19b,id20b)
saveRDS(all.id,'data/allid.rds')
saveRDS(all.idb,'data/allidb.rds')

# there there are lines of code that generate warnings, but they aren't errors. 
# This one is fine: In if (is.na(plot.list)) plot.list <- get.plot(year = year) :the condition has length > 1 and only the first element will be used
# the warnings about no non-missing arguments may require further examination
warnings()

