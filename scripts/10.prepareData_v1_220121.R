# Extract plot data for 2013, 2018, 2019, 2020 and prepare for analysis
rm(list=ls())

source('scripts/11.PW_functions_local-test.R')

## FOR NOW CODE TO EXPORT DATA WITH BRANCHES ALL COMMENTED OUT

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
id13b[id13b$Num==1671,]

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
saveRDS(all.id,'data/allid.Rdata')
saveRDS(all.idb,'data/allidb.Rdata')

# there are lines of code that generate warnings, but they aren't errors. 
# This one is fine: In if (is.na(plot.list)) plot.list <- get.plot(year = year) :the condition has length > 1 and only the first element will be used
# the warnings about no non-missing arguments may require further examination
warnings()

