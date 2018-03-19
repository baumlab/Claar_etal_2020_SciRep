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


phy97.f.c.sediment.VH <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$Dist=="VeryHigh")
phy97.f.c.sediment.M <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$Dist=="HighMed")

phy97.f.c.water.VH <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$Dist=="VeryHigh")
phy97.f.c.water.M <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$Dist=="HighMed")

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




jpeg(filename="figures/pcoa_by_dist_time.jpg", 
     width = 16, height = 8, units="in",res = 300)
par(mfrow=c(2,4))

# Sediment
plot(sediment.VH.bd, hull=F, label=F, 
     main="Sediment-VH", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.VH.bd, sample_data(phy97.f.c.sediment.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.sediment)$field_season)),
       col=timecols)

plot(sediment.M.bd, hull=F, label=F, 
     main="Sediment-M", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.M.bd, sample_data(phy97.f.c.sediment.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.sediment)$field_season)),
       col=timecols)


# Water
plot(water.VH.bd, hull=F, label=F, 
     main="water-VH", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.VH.bd, sample_data(phy97.f.c.water.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.water)$field_season)),
       col=timecols)

plot(water.M.bd, hull=F, label=F, 
     main="water-M", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.M.bd, sample_data(phy97.f.c.water.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.water)$field_season)),
       col=timecols)

boxplot(sediment.VH.bd,col=timecols,ylim=c(0,0.101))
boxplot(water.VH.bd,col=timecols,ylim=c(0,0.101))
boxplot(sediment.M.bd,col=timecols,ylim=c(0,0.101))
boxplot(water.M.bd,col=timecols,ylim=c(0,0.101))


dev.off()

anova(sediment.VH.bd)
anova(sediment.M.bd)
anova(water.VH.bd)
anova(water.M.bd)

phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

phy97.f.c.coral.Peyd.VH <- subset_samples(phy97.f.c.coral.Peyd,sample_data(phy97.f.c.coral.Peyd)$Dist=="VeryHigh")
phy97.f.c.coral.MAeq.VH <- subset_samples(phy97.f.c.coral.MAeq,sample_data(phy97.f.c.coral.MAeq)$Dist=="VeryHigh")
phy97.f.c.coral.Plob.VH <- subset_samples(phy97.f.c.coral.Plob,sample_data(phy97.f.c.coral.Plob)$Dist=="VeryHigh")

phy97.f.c.coral.Peyd.M <- subset_samples(phy97.f.c.coral.Peyd,sample_data(phy97.f.c.coral.Peyd)$Dist=="HighMed")
phy97.f.c.coral.MAeq.M <- subset_samples(phy97.f.c.coral.MAeq,sample_data(phy97.f.c.coral.MAeq)$Dist=="HighMed")
phy97.f.c.coral.Plob.M <- subset_samples(phy97.f.c.coral.Plob,sample_data(phy97.f.c.coral.Plob)$Dist=="HighMed")


Peyd.VH.ufdist <- UniFrac(phy97.f.c.coral.Peyd.VH, weighted=T, normalized=F, parallel=F, fast=T)
Peyd.VH.bd <- betadisper(d=Peyd.VH.ufdist, 
                             group=sample_data(phy97.f.c.coral.Peyd.VH)$field_season,
                             type="centroid", bias.adjust=FALSE)
Plob.VH.ufdist <- UniFrac(phy97.f.c.coral.Plob.VH, weighted=T, normalized=F, parallel=F, fast=T)
Plob.VH.bd <- betadisper(d=Plob.VH.ufdist, 
                         group=sample_data(phy97.f.c.coral.Plob.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Maeq.VH.ufdist <- UniFrac(phy97.f.c.coral.MAeq.VH, weighted=T, normalized=F, parallel=F, fast=T)
Maeq.VH.bd <- betadisper(d=Maeq.VH.ufdist, 
                         group=sample_data(phy97.f.c.coral.MAeq.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)

Peyd.M.ufdist <- UniFrac(phy97.f.c.coral.Peyd.M, weighted=T, normalized=F, parallel=F, fast=T)
Peyd.M.bd <- betadisper(d=Peyd.M.ufdist, 
                         group=sample_data(phy97.f.c.coral.Peyd.M)$field_season,
                         type="centroid", bias.adjust=FALSE)
Plob.M.ufdist <- UniFrac(phy97.f.c.coral.Plob.M, weighted=T, normalized=F, parallel=F, fast=T)
Plob.M.bd <- betadisper(d=Plob.M.ufdist, 
                         group=sample_data(phy97.f.c.coral.Plob.M)$field_season,
                         type="centroid", bias.adjust=FALSE)
Maeq.M.ufdist <- UniFrac(phy97.f.c.coral.MAeq.M, weighted=T, normalized=F, parallel=F, fast=T)
Maeq.M.bd <- betadisper(d=Maeq.M.ufdist, 
                         group=sample_data(phy97.f.c.coral.MAeq.M)$field_season,
                         type="centroid", bias.adjust=FALSE)



plot(Peyd.VH.bd, hull=F, label=F, 
     main="Peyd-VH", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Peyd.VH.bd, sample_data(phy97.f.c.coral.Peyd.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd.VH)$field_season)),
       col=timecols)

plot(Peyd.M.bd, hull=F, label=F, 
     main="Peyd-M", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Peyd.M.bd, sample_data(phy97.f.c.coral.Peyd.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd.VH)$field_season)),
       col=timecols)

plot(Plob.VH.bd, hull=F, label=F, 
     main="Plob-VH", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Plob.VH.bd, sample_data(phy97.f.c.coral.Plob.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob.VH)$field_season)),
       col=timecols)

plot(Plob.M.bd, hull=F, label=F, 
     main="Plob-M", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Plob.M.bd, sample_data(phy97.f.c.coral.Plob.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob.VH)$field_season)),
       col=timecols)

plot(Maeq.VH.bd, hull=F, label=F, 
     main="MAeq-VH", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Maeq.VH.bd, sample_data(phy97.f.c.coral.MAeq.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq.VH)$field_season)),
       col=timecols)

plot(Maeq.M.bd, hull=F, label=F, 
     main="MAeq-M", col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(Maeq.M.bd, sample_data(phy97.f.c.coral.MAeq.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq.VH)$field_season)),
       col=timecols)

boxplot(Maeq.VH.bd,col=timecols,ylim=c(0,0.05))
boxplot(Plob.VH.bd,col=timecols,ylim=c(0,0.05))
boxplot(Peyd.VH.bd,col=timecols,ylim=c(0,0.05))
boxplot(Maeq.M.bd,col=timecols,ylim=c(0,0.05))
boxplot(Plob.M.bd,col=timecols,ylim=c(0,0.05))
boxplot(Peyd.M.bd,col=timecols,ylim=c(0,0.05))
