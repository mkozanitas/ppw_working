# compare seju counts to ladder fuels
rm(list=ls())

# read in seju data summarized by plot
seju.data.p <- read.csv('data/seju_data_p.csv')
dim(seju.data.p)

# read in ladder fuel data
ld <- read.csv('data/mean_ladder_fuel_54.csv')
dim(ld)

seju.data.p$ld <- ld$ld[1:50]
head(seju.data.p)

#plot(seju.data.p$ld~seju.data.p$Total.Number)
plot(ld~Total.Number,data=seju.data.p)
fit <- lm(ld~Total.Number,data=seju.data.p)
abline(fit)
