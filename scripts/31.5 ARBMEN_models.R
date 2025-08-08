library(brms)

set.seed(321)

#dd <- readRDS("brm.ARBMEN.dd.rds")

ddsub <- dd[c(),]
for(i in 1:length(unique(dd$TreeNum))){
  ddtemp <- subset(dd, TreeNum==unique(dd$TreeNum)[i])
  ddtemp <- ddtemp[sample(1:nrow(ddtemp), size=1),]
  ddsub <- rbind(ddsub, ddtemp)
}

multifit_nosubset <- brm(fate3.18 ~ s(d10.17, k=3, by=fac.fsCat) + fac.fsCat + (1|Plot), data=dd,
                         family="categorical",
                         chains = 2,
                         cores = 2,
                         seed=726,
                         iter=50000,
                         refresh=100,
                         control=list(adapt_delta=0.8))

saveRDS(summary(warnings()),paste(local.dir,'/brm.',spName,'.ALL.',uhn,'.MN.Splk3.fate3.18.i',iter,'.WARNINGS.rds',sep=''))
saveRDS(multifit_nosubset,paste(local.dir,'/brm.',spName,'.ALL.',uhn,'.MN.Splk3.fate3.18.i',iter,'.rds',sep=''))
print(summary(multifit_nosubset))

visualizeMultifitBayes(multifit_nosubset,sp='MN.Splk3') 
conditional_effects(multifit_nosubset, categorical=TRUE)


multifit_subset <- brm(fate3.18 ~ s(d10.17, k=3, by=fac.fsCat) + fac.fsCat + (1|Plot), data=ddsub,
                       family="categorical",
                       chains = 2,
                       cores = 2,
                       seed=726,
                       iter=50000,
                       refresh=100,
                       control=list(adapt_delta=0.8))

saveRDS(summary(warnings()),paste(local.dir,'/brm.',spName,'.SUB.',uhn,'.MN.Splk3.fate3.18.i',iter,'.WARNINGS.rds',sep=''))
saveRDS(multifit_subset,paste(local.dir,'/brm.',spName,'.SUB.',uhn,'.MN.Splk3.fate3.18.i',iter,'.rds',sep=''))
print(summary(multifit_subset))

visualizeMultifitBayes(multifit_subset,sp='MN.Splk3') 

multifit_nosubset
multifit_subset

conditional_effects(multifit_nosubset, categorical=TRUE)
conditional_effects(multifit_subset, categorical=TRUE)

# Results look quite similar.