# sandbox

# code based on first round of structure based ordination
# find a couple plots to examine
which.max(ss[,2])
ss[108,]
ps[54*c(1:3),]

table(tAll$Plot.17)
p54 <- which(tAll$Plot.17=='PPW1854' | tAll$Plot.18=='PPW1854' | tAll$Plot.19=='PPW1854' | tAll$Plot.20=='PPW1854')

p54 <- tAll[p54,]
dim(p54)
table(p54$Live.17,useNA = 'always')
table(p54$Live.18,useNA = 'always')
table(p54$Live.19,useNA = 'always')
table(p54$Live.20,useNA = 'always')

table(p54$Topkill.18,p54$Topkill.19)
table(p54$gCrown.17,p54$gCrown.18,p54$gCrown.19)

lseq <- paste(p54$Live.17,p54$Live.18,p54$Live.19,sep='')
gseq <- paste(p54$gCrown.17,p54$gCrown.18,p54$gCrown.19,sep='')
tseq <- paste(p54$Topkill.17,p54$Topkill.18,p54$Topkill.19,sep='')
table(gseq,tseq)

getwd()
write.csv(p54,'data/p54.csv')
ps[54,]

summary(ord)

table(tAll$Plot.17,useNA = 'always')
table(tAll$Plot.18,useNA = 'always')
table(tAll$Plot.19,useNA = 'always')
table(tAll$Plot.20,useNA = 'always')
