options(scipen=999)
library(brms)

# doug fir survival model.
# will probably take a few minutes.
# you might get warnings about "treedepth". this is fine.
# treedepth issues are about computational efficiency, not validity of results.
#d <- read.csv("psemen.csv")
d <- tdat[[1]]

d$fsCat2 <- factor(as.character(d$fsCat))
d$Plot <- factor(as.character(d$Plot))
d$iNum.13 <- factor(as.character(d$iNum.13))

dougfir_splines <- brm(gCrown.18 ~ 
                         s(d10.17, by=fsCat2, k=20)+
                         fsCat2+
                         (1|Plot)+
                         (1|iNum.13),
                       data=d,
                       family="Bernoulli",
                       chains = 2, cores = 2, seed=237, 
                       #backend="cmdstanr",
                       control=list(adapt_delta=0.99))

# quick viz
conditional_effects(dougfir_splines)

# ARCMAN
d <- tdat[[2]]

d$fsCat2 <- factor(as.character(d$fsCat))
d$Plot <- factor(as.character(d$Plot))
d$iNum.13 <- factor(as.character(d$iNum.13))

manzanita_splines <- brm(gCrown.18 ~ 
                         s(d10.17, by=fsCat2, k=20)+
                         fsCat2+
                         (1|Plot)+
                         (1|iNum.13),
                       data=d,
                       family="Bernoulli",
                       chains = 2, cores = 2, seed=237, 
                       #backend="cmdstanr",
                       control=list(adapt_delta=0.99))

# quick viz
conditional_effects(manzanita_splines)

# multinomial model.
# will probably take a few hours,
# or you can just download my attached .RDS file and run:
multifit1 <- readRDS("data/multifit1.RDS")
d2 <- read.csv("ehro.csv")

d2$fsCat2 <- factor(as.character(d2$fsCat))
d2$Plot <- factor(as.character(d2$Plot))
d2$iNum.13 <- factor(as.character(d2$iNum.13))
d2$fatefac <- factor(as.character(d2$fate3.18))

multifit1 <- brm(fatefac ~ s(d10.17, k=20, by=interaction(fsCat2, Species)) + fsCat2*Species + (1|Plot) + (1|iNum.13), data=d2,
                 family="categorical", chains = 2, cores = 2, seed=726, 
                 #backend="cmdstanr",
                 refresh=100,
                 control=list(adapt_delta=0.95))
summary(multifit1)

# canned conditional_effects() isn't great for multinomial models
# the bayesplot package offers lots of options with great tutorials: https://mc-stan.org/bayesplot/

