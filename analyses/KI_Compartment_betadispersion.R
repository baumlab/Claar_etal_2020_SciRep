# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")

phy97.f.c.sediment <- subset_samples(phy97.f.c.sediment,field_season!="KI2015c")

col <- c(KI2014 = "#2c7fb8", KI2015a_Pre = "#7fcdbb", KI2015a_Post = "#253494", KI2015b = "#41b6c4")


sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd <- betadisper(d=sediment.ufdist, group=sample_data(phy97.f.c.sediment)$field_season,
                          type="centroid", bias.adjust=FALSE)

water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, normalized=F, parallel=F, fast=T)
water.bd <- betadisper(d=water.ufdist, group=sample_data(phy97.f.c.water)$field_season,
                       type="centroid", bias.adjust=FALSE)


jpeg(filename="figures/sediment_water_betadisper.jpg", 
     width = 8, height = 4, units="in",res = 300)
par(mfrow=c(1,2))
plot(sediment.bd, hull=F, label=F, 
     main="Sediment", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd, sample_data(phy97.f.c.sediment)$field_season,  draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment)$field_season),
       col=col)

plot(water.bd, hull=F, label=F, 
     main="Water", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd, sample_data(phy97.f.c.water)$field_season,  draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.water)$field_season),
       col=col)


dev.off()

anova(sediment.bd)
anova(water.bd)

jpeg(filename="figures/sediment_water_betadisper_box.jpg", width = 9, height = 4, units="in",res = 300)
par(mfrow=c(1,2),mar=c(3,4,0.5,0.5))
boxplot(sediment.bd,col=col)
boxplot(water.bd,col=col)
dev.off()
