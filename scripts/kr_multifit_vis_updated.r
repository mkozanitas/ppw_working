library(ggplot2)
library(patchwork)

# set up grid for conditional mean predictions
d10.17_grid <- seq(from=min(dd$d10.17), to=max(dd$d10.17), length.out=300)
fac.fsCat_grid <- unique(dd$fac.fsCat)
predgrid <- expand.grid(d10.17_grid, fac.fsCat_grid)
names(predgrid) <- c("d10.17", "fac.fsCat")

# get posterior samples of conditional mean, ignoring random effects
pep <- posterior_epred(object=multifit1, newdata=predgrid, re_formula=NA)

# summarize posterior samples
pepmeans <- as.data.frame(apply(pep, c(2,3), mean))
names(pepmeans) <- paste(names(pepmeans), "_mean", sep="")

pep2.5 <- as.data.frame(apply(pep, c(2,3), quantile, 0.025))
names(pep2.5) <- paste(names(pep2.5), "_q2.5", sep="")

pep97.5 <- as.data.frame(apply(pep, c(2,3), quantile, 0.975))
names(pep97.5) <- paste(names(pep97.5), "_q97.5", sep="")

predgrid <- cbind(predgrid, pepmeans, pep2.5, pep97.5)
predgrid$fac.fsCat <- factor(predgrid$fac.fsCat, levels=c("0.U", "12.LM", "3.H"))

mortplot <- ggplot(data=predgrid, aes(x=d10.17, y=DN_mean, group=fac.fsCat, color=fac.fsCat))+
  ggtitle(spName)+
  geom_line()+
  geom_ribbon(aes(ymin=DN_q2.5, ymax=DN_q97.5, group=fac.fsCat, fill=fac.fsCat, color=NULL), alpha=0.2)+
  scale_color_viridis_d(name="Fire\nSeverity",
                        labels=c("Unburned", "Low/Medium", "High"))+
  scale_fill_viridis_d(name="Fire\nSeverity",
                        labels=c("Unburned", "Low/Medium", "High"))+
  xlab("Log(Diameter at Breast Height [cm])")+ ##### units correct?
  ylab("P(Mortality)")+
  theme_bw()+
  theme(legend.position="none",
        aspect.ratio=1)

resprplot <- ggplot(data=predgrid, aes(x=d10.17, y=DR_mean, group=fac.fsCat, color=fac.fsCat))+
  geom_line()+
  geom_ribbon(aes(ymin=DR_q2.5, ymax=DR_q97.5, group=fac.fsCat, fill=fac.fsCat, color=NULL), alpha=0.2)+
  scale_color_viridis_d(name="Fire\nSeverity",
                        labels=c("Unburned", "Low/Medium", "High"))+
  scale_fill_viridis_d(name="Fire\nSeverity",
                       labels=c("Unburned", "Low/Medium", "High"))+
  xlab("Log(Diameter at Breast Height [cm])")+ ##### units correct?
  ylab("P(Resprout)")+
  theme_bw()+
  theme(legend.position="none",
        aspect.ratio=1)

gcplot <- ggplot(data=predgrid, aes(x=d10.17, y=GC_mean, group=fac.fsCat, color=fac.fsCat))+
  geom_line()+
  geom_ribbon(aes(ymin=GC_q2.5, ymax=GC_q97.5, group=fac.fsCat, fill=fac.fsCat, color=NULL), alpha=0.2)+
  scale_color_viridis_d(name="Fire\nSeverity",
                        labels=c("Unburned", "Low/Medium", "High"))+
  scale_fill_viridis_d(name="Fire\nSeverity",
                       labels=c("Unburned", "Low/Medium", "High"))+
  xlab("Log(Diameter at Breast Height [cm])")+ ##### units correct?
  ylab("P(Green Crown)")+
  theme_bw()+
  theme()

mortplot+resprplot+gcplot




