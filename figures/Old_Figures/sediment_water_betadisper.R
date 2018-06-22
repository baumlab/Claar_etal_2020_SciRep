# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)

sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd <- betadisper(d=sediment.ufdist, group=sample_data(phy97.f.c.sediment)$field_season,
                          type="centroid", bias.adjust=FALSE)

water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, normalized=F, parallel=F, fast=T)
water.bd <- betadisper(d=water.ufdist, group=sample_data(phy97.f.c.water)$field_season,
                       type="centroid", bias.adjust=FALSE)

phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

Peyd.ufdist <- UniFrac(phy97.f.c.coral.Peyd, weighted=T, normalized=F, parallel=F, fast=T)
Peyd.bd <- betadisper(d=Peyd.ufdist, group=sample_data(phy97.f.c.coral.Peyd)$field_season,
                      type="centroid", bias.adjust=FALSE)
Peyd.bd.site <- betadisper(d=Peyd.ufdist, group=sample_data(phy97.f.c.coral.Peyd)$site,
                           type="centroid", bias.adjust=FALSE)

MAeq.ufdist <- UniFrac(phy97.f.c.coral.MAeq, weighted=T, normalized=F, parallel=F, fast=T)
MAeq.bd <- betadisper(d=MAeq.ufdist, group=sample_data(phy97.f.c.coral.MAeq)$field_season,
                      type="centroid", bias.adjust=FALSE)
MAeq.bd.site <- betadisper(d=MAeq.ufdist, group=sample_data(phy97.f.c.coral.MAeq)$site,
                           type="centroid", bias.adjust=FALSE)

Plob.ufdist <- UniFrac(phy97.f.c.coral.Plob, weighted=T, normalized=F, parallel=F, fast=T)
Plob.bd <- betadisper(d=Plob.ufdist, group=sample_data(phy97.f.c.coral.Plob)$field_season,
                      type="centroid", bias.adjust=FALSE)
Plob.bd.site <- betadisper(d=Plob.ufdist, group=sample_data(phy97.f.c.coral.Plob)$site,
                           type="centroid", bias.adjust=FALSE)

jpeg(filename="figures/sediment_water_betadisper.jpg", 
     width = 20, height = 8, units="in",res = 300)
par(mfrow=c(2,5))
plot(sediment.bd, hull=F, label=F, 
    main="Sediment", 
    # main=NA,
     col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd, sample_data(phy97.f.c.sediment)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=c("Before","Post-Storm","After"),
       #legend=levels(sample_data(phy97.f.c.sediment)$field_season),
       col=timecols)

plot(water.bd, hull=F, label=F, 
    main="Water",
    # main=NA,
     col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd, sample_data(phy97.f.c.water)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=c("Before","Post-Storm","After"),
       #legend=levels(sample_data(phy97.f.c.water)$field_season),
       col=timecols)

plot(Peyd.bd, hull=F, label=F, 
     main="Pocillopora eydouxi",
     # main=NA,
     col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Peyd.bd, sample_data(phy97.f.c.coral.Peyd)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=c("Before","Pre-Storm","Post-Storm","After"),
       #legend=levels(sample_data(phy97.f.c.water)$field_season),
       col=timecols)

plot(MAeq.bd, hull=F, label=F, 
     main="Montipora aequituburculata",
     # main=NA,
     col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(MAeq.bd, sample_data(phy97.f.c.coral.MAeq)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=c("Before","Pre-Storm","Post-Storm","After"),
       #legend=levels(sample_data(phy97.f.c.water)$field_season),
       col=timecols)


plot(Plob.bd, hull=F, label=F, 
     main="Porites lobata",
     # main=NA,
     col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Plob.bd, sample_data(phy97.f.c.coral.Plob)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=c("Before","Pre-Storm","Post-Storm","After"),
       #legend=levels(sample_data(phy97.f.c.water)$field_season),
       col=timecols)


boxplot(sediment.bd,col=timecols,names=c("Before","Post-Storm","After"))
boxplot(water.bd,col=timecols,names=c("Before","Post-Storm","After"))
boxplot(Peyd.bd,col=timecols,names=c("Before","Pre-Storm","Post-Storm","After"))
boxplot(MAeq.bd,col=timecols,names=c("Before","Pre-Storm","Post-Storm","After"))
boxplot(Plob.bd,col=timecols,names=c("Before","Pre-Storm","Post-Storm","After"))


dev.off()

anova(sediment.bd)
anova(water.bd)
anova(Peyd.bd)
anova(MAeq.bd)
anova(Plob.bd)
