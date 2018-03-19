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

sediment.bd.site <- betadisper(d=sediment.ufdist, group=sample_data(phy97.f.c.sediment)$site,
                               type="centroid", bias.adjust=FALSE)
water.bd.site <- betadisper(d=water.ufdist, group=sample_data(phy97.f.c.water)$site,
                            type="centroid", bias.adjust=FALSE)


phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

Peyd.ufdist <- UniFrac(phy97.f.c.coral.Peyd, weighted=T, normalized=F, parallel=F, fast=T)
Peyd.bd <- betadisper(d=Peyd.ufdist, group=sample_data(phy97.f.c.coral.Peyd)$field_season,
                          type="centroid", bias.adjust=FALSE)
Peyd.bd.site <- betadisper(d=Peyd.ufdist, group=sample_data(phy97.f.c.coral.Peyd)$site,
                               type="centroid", bias.adjust=FALSE)
Peyd.bd.dist <- betadisper(d=Peyd.ufdist, group=sample_data(phy97.f.c.coral.Peyd)$Dist,
                           type="centroid", bias.adjust=FALSE)

MAeq.ufdist <- UniFrac(phy97.f.c.coral.MAeq, weighted=T, normalized=F, parallel=F, fast=T)
MAeq.bd <- betadisper(d=MAeq.ufdist, group=sample_data(phy97.f.c.coral.MAeq)$field_season,
                      type="centroid", bias.adjust=FALSE)
MAeq.bd.site <- betadisper(d=MAeq.ufdist, group=sample_data(phy97.f.c.coral.MAeq)$site,
                           type="centroid", bias.adjust=FALSE)
MAeq.bd.dist <- betadisper(d=MAeq.ufdist, group=sample_data(phy97.f.c.coral.MAeq)$Dist,
                           type="centroid", bias.adjust=FALSE)

Plob.ufdist <- UniFrac(phy97.f.c.coral.Plob, weighted=T, normalized=F, parallel=F, fast=T)
Plob.bd <- betadisper(d=Plob.ufdist, group=sample_data(phy97.f.c.coral.Plob)$field_season,
                      type="centroid", bias.adjust=FALSE)
Plob.bd.site <- betadisper(d=Plob.ufdist, group=sample_data(phy97.f.c.coral.Plob)$site,
                           type="centroid", bias.adjust=FALSE)
Plob.bd.dist <- betadisper(d=Plob.ufdist, group=sample_data(phy97.f.c.coral.Plob)$Dist,
                           type="centroid", bias.adjust=FALSE)

coral.ufdist <- UniFrac(phy97.f.c.coral, weighted=T, normalized=F, parallel=F, fast=T)
coral.bd <- betadisper(d=coral.ufdist, group=sample_data(phy97.f.c.coral)$field_season,
                      type="centroid", bias.adjust=FALSE)
coral.bd.site <- betadisper(d=coral.ufdist, group=sample_data(phy97.f.c.coral)$site,
                           type="centroid", bias.adjust=FALSE)
coral.bd.dist <- betadisper(d=coral.ufdist, group=sample_data(phy97.f.c.coral)$Dist,
                           type="centroid", bias.adjust=FALSE)
coral.bd.species <- betadisper(d=coral.ufdist, group=sample_data(phy97.f.c.coral)$Coral_Species,
                            type="centroid", bias.adjust=FALSE)


jpeg(filename="figures/pcoa_by_site.jpg", 
     width = 20, height = 8, units="in",res = 300)
par(mfrow=c(2,5))

# Sediment
plot(sediment.bd.site, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.site, sample_data(phy97.f.c.sediment)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment)$site),
       col=sitecols)


# Water
plot(water.bd.site, hull=F, label=F, 
     main="Water", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd.site, sample_data(phy97.f.c.water)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.water)$site),
       col=sitecols)

# Pocillopora
plot(Peyd.bd.site, hull=F, label=F, 
     main="Pocillopora eydouxi", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Peyd.bd.site, sample_data(phy97.f.c.coral.Peyd)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.coral.Peyd)$site),
       col=sitecols)


# Montipora
plot(MAeq.bd.site, hull=F, label=F, 
     main="Montipora aequituburculata", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(MAeq.bd.site, sample_data(phy97.f.c.coral.MAeq)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.coral.MAeq)$site),
       col=sitecols)

# Porites
plot(Plob.bd.site, hull=F, label=F, 
     main="Porites lobata", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Plob.bd.site, sample_data(phy97.f.c.coral.Plob)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.coral.Plob)$site),
       col=sitecols)

boxplot(sediment.bd.site,col=sitecols,names=c("Before","Post-Storm","After"))
boxplot(water.bd.site,col=sitecols,names=c("Before","Post-Storm","After"))
boxplot(Peyd.bd.site,col=sitecols,names=c("Before","Pre-Storm","Post-Storm","After"))
boxplot(MAeq.bd.site,col=sitecols,names=c("Before","Pre-Storm","Post-Storm","After"))
boxplot(Plob.bd.site,col=sitecols,names=c("Before","Pre-Storm","Post-Storm","After"))


dev.off()

anova(sediment.bd.site)
anova(water.bd.site)
anova(Peyd.bd.site)
anova(MAeq.bd.site)
anova(Plob.bd.site)


anova(Peyd.bd.dist)
anova(MAeq.bd.dist)
anova(Plob.bd.dist)


plot(coral.bd.species, hull=F, label=T, 
     main="coral",
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(coral.bd.species, sample_data(phy97.f.c.coral)$Coral_Species,  draw = c("polygon"),
         alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(sample_data(phy97.f.c.coral)$Coral_Species))

