## read in niche values and subset to PWD species
cn <- read.csv('/Users/david/My Drive/My_Drive_Cloud/Drive-Projects/CalDiversity/ClimaticNicheAnalysis/CA_CNS_2017_draft_offical/climate_optima.csv')

sp <- read.csv('data/pwd-species-codes.csv')

cn_pwd <- cn[which(cn$spp %in% sp$GS),]
cn_pwd$spcode <- sp$SpCode[match(cn_pwd$spp,sp$GS)]
dim(cn_pwd)
head(cn_pwd)

write.csv(cn_pwd,'data/pwd-species-niche.csv',row.names = FALSE)
