#Summarizing Basal Area using 2018 data and filling in missing trees using 

rm(list=ls())
tAll <- read.csv('data/tAll.csv')

library(dplyr)
library(tidyverse)

AbSp <- c('ARBMEN','ARCMAN','HETARB','PSEMEN','QUEAGR','QUEDOU','QUEGAR','QUEKEL','UMBCAL')

tAll$BAmerge <-ifelse(is.na(tAll$Basal.Area.18),tAll$Basal.Area.13,tAll$Basal.Area.18)
tAll$Sp.Type <- paste(tAll$Type.13,tAll$Species.13,sep='_')

tAll %>%
  filter(Type.13=='TR') %>% 
  group_by(Species.13,fate.18) %>%
  summarize(BAspecies=sum(BAmerge,na.rm = T)) %>% 
  pivot_wider(names_from = fate.18,values_from = BAspecies) %>% 
  mutate(BAsum=sum(DN,DR,LN,LR,na.rm=T),BAdelta=-1*sum(DN,DR,na.rm=T)) %>% 
  filter(Species.13 %in% AbSp) 

tAll %>%
  group_by(Sp.Type,fate.18) %>%
  summarize(BAspecies=sum(BAmerge,na.rm = T)) %>% 
  pivot_wider(names_from = fate.18,values_from = BAspecies) %>% 
  mutate(BAsum=sum(DN,DR,LN,LR,na.rm=T),BAdelta=-1*sum(DN,DR,na.rm=T))

table(tAll$Type.13,tAll$Type.18,useNA = 'always')
table(tAll$Type.18,tAll$Type.19,useNA = 'always')
table(tAll$Type.19,tAll$Type.20,useNA = 'always')



tAll %>%
  filter(Type.18=='TR') %>% 
  group_by(fate.18) %>%
  summarize(BAspecies=sum(BAmerge,na.rm = T)) %>% 
  pivot_wider(names_from = fate.18,values_from = BAspecies) %>% 
  mutate(BAsum=sum(DN,DR,LN,LR,na.rm=T),BAdelta=-1*sum(DN,DR,na.rm=T))
