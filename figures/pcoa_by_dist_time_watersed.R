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

sediment.VH.ufdist <- UniFrac(phy97.f.c.sediment.VH, weighted=T, normalized=F, parallel=F, fast=T)
water.VH.ufdist <- UniFrac(phy97.f.c.water.VH, weighted=T, normalized=F, parallel=F, fast=T)
sediment.VH.bd <- betadisper(d=sediment.VH.ufdist, 
                             group=sample_data(phy97.f.c.sediment.VH)$field_season,
                             type="centroid", bias.adjust=FALSE)
water.VH.bd <- betadisper(d=water.VH.ufdist, group=sample_data(phy97.f.c.water.VH)$field_season,
                          type="centroid", bias.adjust=FALSE)

sediment.M.ufdist <- UniFrac(phy97.f.c.sediment.M, weighted=T, normalized=F, parallel=F, fast=T)
water.M.ufdist <- UniFrac(phy97.f.c.water.M, weighted=T, normalized=F, parallel=F, fast=T)
sediment.M.bd <- betadisper(d=sediment.M.ufdist, 
                            group=sample_data(phy97.f.c.sediment.M)$field_season,
                            type="centroid", bias.adjust=FALSE)
water.M.bd <- betadisper(d=water.M.ufdist, group=sample_data(phy97.f.c.water.M)$field_season,
                         type="centroid", bias.adjust=FALSE)

jpeg(filename="figures/pcoa_by_dist_time_watersed.jpg", 
     width = 8, height = 8, units="in",res = 300)
par(mfrow=c(2,2))

# Sediment
plot(sediment.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title("Sediment", line = 1)
ordihull(sediment.VH.bd, sample_data(phy97.f.c.sediment.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("bottomright", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.sediment)$field_season)),
       col=timecols)
mtext("VERY HIGH DISTURBANCE",side=3,line=3)

plot(sediment.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title("Sediment", line = 1)
ordihull(sediment.M.bd, sample_data(phy97.f.c.sediment.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.sediment)$field_season)),
       col=timecols)
mtext("MEDIUM DISTURBANCE",side=3,line=3)


# Water
plot(water.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title("Water", line = 1)
ordihull(water.VH.bd, sample_data(phy97.f.c.water.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("bottomright", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.water)$field_season)),
       col=timecols)

plot(water.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title("Water", line = 1)
ordihull(water.M.bd, sample_data(phy97.f.c.water.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.water)$field_season)),
       col=timecols)

dev.off()