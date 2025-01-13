## troubleshooting gCrown data
xx <- id18b[-which(id18b$Type=='TS'),]
table(xx$Topkill,xx$gCrown,useNA='always')
table(xx$Survival,xx$gCrown,useNA='always')
table(xx$Survival,xx$gCrown,xx$Topkill,useNA='always')

# it seems that all Topkill==1 should be gCrown=0
T1gNA <- which(xx$Topkill==1 & is.na(xx$gCrown))
length(T1gNA)
head(xx[T1gNA,])

# what about Topkill==0 and gCrown=NA
T0gNA <- which(xx$Topkill==0 & is.na(xx$gCrown))
length(T0gNA)
head(xx[T0gNA,])
tail(id18b[T0gNA,])
