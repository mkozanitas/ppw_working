## read in niche values and subset to PWD species
cn <- read.csv('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/CalDiversity/ClimaticNicheAnalysis/CA_CNS_2017_draft_offical/climate_optima.csv')
head(cn)

sp <- read.csv('data/all-spp-names.csv')
head(sp)
cn_pwd <- cn[which(cn$spp %in% sp$species),]
cn_pwd$sp_code <- sp$sp_code[match(cn_pwd$spp,sp$species)]
dim(cn_pwd)
head(cn_pwd)

write.csv(cn_pwd,'data/pwd-species-niche.csv',row.names = FALSE)
