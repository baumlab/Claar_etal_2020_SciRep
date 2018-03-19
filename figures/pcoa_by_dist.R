# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)

sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, normalized=F, parallel=F, fast=T)
water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd.dist <- betadisper(d=sediment.ufdist, group=sample_data(phy97.f.c.sediment)$Dist,
                               type="centroid", bias.adjust=FALSE)
water.bd.dist <- betadisper(d=water.ufdist, group=sample_data(phy97.f.c.water)$Dist,
                            type="centroid", bias.adjust=FALSE)


phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

Peyd.ufdist <- UniFrac(phy97.f.c.coral.Peyd, weighted=T, normalized=F, parallel=F, fast=T)
Peyd.bd.dist <- betadisper(d=Peyd.ufdist, group=sample_data(phy97.f.c.coral.Peyd)$Dist,
                           type="centroid", bias.adjust=FALSE)

MAeq.ufdist <- UniFrac(phy97.f.c.coral.MAeq, weighted=T, normalized=F, parallel=F, fast=T)
MAeq.bd.dist <- betadisper(d=MAeq.ufdist, group=sample_data(phy97.f.c.coral.MAeq)$Dist,
                           type="centroid", bias.adjust=FALSE)

Plob.ufdist <- UniFrac(phy97.f.c.coral.Plob, weighted=T, normalized=F, parallel=F, fast=T)
Plob.bd.dist <- betadisper(d=Plob.ufdist, group=sample_data(phy97.f.c.coral.Plob)$Dist,
                           type="centroid", bias.adjust=FALSE)

coral.ufdist <- UniFrac(phy97.f.c.coral, weighted=T, normalized=F, parallel=F, fast=T)
coral.bd.dist <- betadisper(d=coral.ufdist, group=sample_data(phy97.f.c.coral)$Dist,
                            type="centroid", bias.adjust=FALSE)
coral.bd.species <- betadisper(d=coral.ufdist, group=sample_data(phy97.f.c.coral)$Coral_Species,
                               type="centroid", bias.adjust=FALSE)


jpeg(filename="figures/pcoa_by_dist.jpg", 
     width = 20, height = 8, units="in",res = 300)
par(mfrow=c(2,5))

# Sediment
plot(sediment.bd.dist, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.dist, sample_data(phy97.f.c.sediment)$Dist,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.sediment)$Dist)),
       col=sitecols)


# Water
plot(water.bd.dist, hull=F, label=F, 
     main="Water", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd.dist, sample_data(phy97.f.c.water)$Dist,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.water)$Dist)),
       col=sitecols)

# Pocillopora
plot(Peyd.bd.dist, hull=F, label=F, 
     main="Pocillopora eydouxi", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Peyd.bd.dist, sample_data(phy97.f.c.coral.Peyd)$Dist,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd)$Dist)),
       col=sitecols)


# Montipora
plot(MAeq.bd.dist, hull=F, label=F, 
     main="Montipora aequituburculata", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(MAeq.bd.dist, sample_data(phy97.f.c.coral.MAeq)$Dist,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq)$Dist)),
       col=sitecols)

# Porites
plot(Plob.bd.dist, hull=F, label=F, 
     main="Porites lobata", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Plob.bd.dist, sample_data(phy97.f.c.coral.Plob)$Dist,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob)$Dist)),
       col=sitecols)

boxplot(sediment.bd.dist,col=sitecols,ylim=c(0,0.12))
boxplot(water.bd.dist,col=sitecols,ylim=c(0,0.12))
boxplot(Peyd.bd.dist,col=sitecols,ylim=c(0,0.12))
boxplot(MAeq.bd.dist,col=sitecols,ylim=c(0,0.12))
boxplot(Plob.bd.dist,col=sitecols,ylim=c(0,0.12))


dev.off()

anova(sediment.bd.dist)
anova(water.bd.dist)
anova(Peyd.bd.dist)
anova(MAeq.bd.dist)
anova(Plob.bd.dist)

