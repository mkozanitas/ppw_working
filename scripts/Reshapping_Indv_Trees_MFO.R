# Intent: Reshape all.id dataframe made it examineData.R 
# M.F. Oldfather
# Date Created: 20220601
# Date Last Edited: 20220601

library(tidyverse)

# Create all.id dataframe made in examineData.R (fates (not size) of individual trees)
source("scripts/examineData.R")
all.id  # list of 4 dataframes (2013, 2018, 2019, 2020)
str(all.id) # 

# want to combine the list together
do.call("rbind", all.id)
# the columns do not seem match - let's look into that
colnames(all.id[[1]]) # 21
colnames(all.id[[2]]) # 25
colnames(all.id[[3]])
colnames(all.id[[4]])
# the latter list elements (not 2018) have the columns "pattern, DT, TG, DG --> can we add this columns into 2018 as we only tagged certain types of individuals? Let's just guess and say they are all 1 for those columns in 2018

# create data frame to play with
indivs <- all.id

# add in missing 2018 columns with all ones 
indivs[[1]]$pattern <- "11111111"
indivs[[1]]$DT <- 1
indivs[[1]]$TG <- 1
indivs[[1]]$DG <- 1

# rearrange so columns are in the in same order as later years
indivs[[1]] <- indivs[[1]][c(1:17,22:25, 18:21)]
indivs

# merge list 'years' together
indivs_all_years <- do.call("rbind", indivs)
indivs_all_years

# reshape data so the columns are years
indivs_all_years %>% 
  pivot_wider(id_cols = Num, names_from = Year, values_from = "pattern")
