library(betareg)
library(statmod)
op <- readRDS(file.choose())
head(op)
dim(op)

op$tav <- (op$tmax+op$tmin)/2
opscore <- (op$ba_oak-op$ba_pine)/(op$ba_oak+op$ba_pine)
opscaled <- 0.5+opscore/2
plot(opscore~tmax,data=op)
plot(opscore~tmin,data=op)
plot(opscore~tmax,data=op)
plot(opscore~tav,data=op)

fit <- betareg(opscaled ~ tmax, data = op)
summary(fit)

op_test <- op
op_test$pop <- predict(fit,newdata=op_test)
plot(opscaled~tmax,data=op)
points(pop~tmax,data=op_test,pch=19,col='red')



fpscore <- (op$ba_fagaceae-op$ba_pinaceae)/(op$ba_fagaceae+op$ba_pinaceae)
fpscaled <- 0.5+fpscore/2
plot(fpscore~tmax,data=op)
plot(fpscore~tmin,data=op)


