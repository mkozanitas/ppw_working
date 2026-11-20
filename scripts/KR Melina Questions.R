KR cleanup - questions for Melina

Hectares:
  Is PercLivingCanopy the same as Canopy.Percent in plot data?
  Does Senesced have an equivalent in plot data?
  
  KR10.5 - line 391
## MELINA - should we remove outliers or move them to nearest subplot?

L27 # Column named X is a note in wrong column PPW1340 9758 "changed from umbcal due to 2018 fieldnote- mkqc2025"


Plots:
  KR26_20
L58 # 9/24/24 - Melina: There are 955 NA for bSprout in 2013. Is that okay? Make them all 0?

L193
# This individual is a sapling with no size data in 2013 so it can't be used for analysis. It can be left in the data as it will not causes any problems at the analysis step, since it doesn't have a size. And since it was killed I would trust the 2013 species ID more than 2018. Changing the 2018 id here to match 2013
## EDIT AFTER MELINA FIX
all.id[[2]]$Species[which(all.id[[1]]$Num==4072)] <- 'QUEAGR'

L422   # NEED TO REVIEW ALL THIS WITH MELINA

L459   # output suspicious transitions- can switch Type and year here to look at other combos
# MELINA_CHECK

L665
# a few errors in Types - fix these in csvs?
## EDIT AFTER MELINA FIX

L695
# assign cumulative tag pulled variable - run again after fixing 2020 - Melina
tAllh$Tag.Pulled.Cum <- 0
tAllh$Tag.Pulled.Cum[which(tAllh$Tag.Pulled.18==1)] <- 1
tAllh$Tag.Pulled.Cum[which(tAllh$Tag.Pulled.19==1)] <- 1
tAllh$Tag.Pulled.Cum[which(tAllh$Tag.Pulled.20=="1")] <- 1
table(tAllh$Tag.Pulled.Cum)

## EDIT AFTER MELINA FIX
tp20comments <- which(tAllh$Tag.Pulled.20 %in% c('LOST','NEW SA','NOT FOUND','NOT FOUND IN 19 OR 20','NOT FOUND, DOSER KILL','TAG MOVED TO RSPRT','TYPE IN 19, NOT TR'))
write.csv(tAllh[tp20comments,c('Plot','Num','Tag.Pulled.20')],'data/tp20problems.csv')

L706
# Construct a cumulative Quad, X, Y for plot stems, using most recent values
## NEED TO VET THIS WITH MELINA

L733
# CALCULATE UTM for all plot indvs - TO DO

