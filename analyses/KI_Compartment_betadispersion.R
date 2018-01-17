# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")

phy97.f.c.sediment <- subset_samples(phy97.f.c.sediment,field_season!="KI2015c")

col <- c(KI2014 = "#2c7fb8", KI2015a_Pre = "#7fcdbb", KI2015a_Post = "#253494", KI2015b = "#41b6c4")
sitecols <- c("8" = "#1b9e77","30"="#d95f02","35" ="#7570b3", "27"="#66a61e")

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

sediment.bd.site <- betadisper(d=sediment.ufdist, group=sample_data(phy97.f.c.sediment)$site,
                          type="centroid", bias.adjust=FALSE)

plot(sediment.bd.site, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.site, sample_data(phy97.f.c.sediment)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment)$site),
       col=sitecols)
anova(sediment.bd.site)


water.bd.site <- betadisper(d=water.ufdist, group=sample_data(phy97.f.c.water)$site,
                               type="centroid", bias.adjust=FALSE)

plot(water.bd.site, hull=F, label=F, 
     main="water", col=sitecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd.site, sample_data(phy97.f.c.water)$site,  draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.water)$site),
       col=sitecols)
anova(water.bd.site)

jpeg(filename="figures/sediment_betadisper_bysite.jpg", width = 18, height = 4, units="in",res = 300)
par(mfrow=c(1,4))
phy97.f.c.sediment.site27 <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$site==27, prune=TRUE)
sediment.ufdist.site27 <- UniFrac(phy97.f.c.sediment.site27, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd.site27 <- betadisper(d=sediment.ufdist.site27, 
                                 group=sample_data(phy97.f.c.sediment.site27)$field_season,
                                 type="centroid", bias.adjust=FALSE)
plot(sediment.bd.site27, hull=F, label=F, 
     main="Sediment - 27", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.site27, sample_data(phy97.f.c.sediment.site27)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment.site27)$field_season),
       col=col)
anova(sediment.bd.site27)

phy97.f.c.sediment.site30 <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$site==30, prune=TRUE)
sediment.ufdist.site30 <- UniFrac(phy97.f.c.sediment.site30, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd.site30 <- betadisper(d=sediment.ufdist.site30, 
                                 group=sample_data(phy97.f.c.sediment.site30)$field_season,
                                 type="centroid", bias.adjust=FALSE)
plot(sediment.bd.site30, hull=F, label=F, 
     main="Sediment - 30", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.site30, sample_data(phy97.f.c.sediment.site30)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment.site30)$field_season),
       col=col)
anova(sediment.bd.site30)

phy97.f.c.sediment.site35 <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$site==35, prune=TRUE)
sediment.ufdist.site35 <- UniFrac(phy97.f.c.sediment.site35, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd.site35 <- betadisper(d=sediment.ufdist.site35, 
                                 group=sample_data(phy97.f.c.sediment.site35)$field_season,
                                 type="centroid", bias.adjust=FALSE)
plot(sediment.bd.site35, hull=F, label=F, 
     main="Sediment - 35", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.site35, sample_data(phy97.f.c.sediment.site35)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment.site35)$field_season),
       col=col)
anova(sediment.bd.site35)

phy97.f.c.sediment.site8 <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$site==8, prune=TRUE)
sediment.ufdist.site8 <- UniFrac(phy97.f.c.sediment.site8, weighted=T, normalized=F, parallel=F, fast=T)
sediment.bd.site8 <- betadisper(d=sediment.ufdist.site8, 
                                 group=sample_data(phy97.f.c.sediment.site8)$field_season,
                                 type="centroid", bias.adjust=FALSE)
plot(sediment.bd.site8, hull=F, label=F, 
     main="Sediment - 8", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(sediment.bd.site8, sample_data(phy97.f.c.sediment.site8)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.sediment.site8)$field_season),
       col=col)
anova(sediment.bd.site8)
dev.off()

jpeg(filename="figures/sediment_betadisper_boxplot_bysite.jpg", width = 18, height = 4, units="in",res = 300)
par(mfrow=c(1,4))
boxplot(sediment.bd.site27,col=col,main="Site 27")
boxplot(sediment.bd.site30,col=col,main="Site 30")
boxplot(sediment.bd.site35,col=col,main="Site 35")
boxplot(sediment.bd.site8,col=col,main="Site 8")
dev.off()




jpeg(filename="figures/water_betadisper_bysite.jpg", width = 13.5, height = 4, units="in",res = 300)
par(mfrow=c(1,3))
# phy97.f.c.water.site27 <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$site==27, prune=TRUE)
# water.ufdist.site27 <- UniFrac(phy97.f.c.water.site27, weighted=T, normalized=F, parallel=F, fast=T)
# water.bd.site27 <- betadisper(d=water.ufdist.site27, 
#                                  group=sample_data(phy97.f.c.water.site27)$field_season,
#                                  type="centroid", bias.adjust=FALSE)
# plot(water.bd.site27, hull=F, label=F, 
#      main="water - 27", col=col,
#      xlab="PCoA 1", ylab="PCoA 2", sub="")
# ordihull(water.bd.site27, sample_data(phy97.f.c.water.site27)$field_season,  
#          draw = c("polygon"),
#          col = col, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(sample_data(phy97.f.c.water.site27)$field_season),
#        col=col)
# anova(water.bd.site27)

phy97.f.c.water.site30 <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$site==30, prune=TRUE)
water.ufdist.site30 <- UniFrac(phy97.f.c.water.site30, weighted=T, normalized=F, parallel=F, fast=T)
water.bd.site30 <- betadisper(d=water.ufdist.site30, 
                                 group=sample_data(phy97.f.c.water.site30)$field_season,
                                 type="centroid", bias.adjust=FALSE)
plot(water.bd.site30, hull=F, label=F, 
     main="water - 30", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd.site30, sample_data(phy97.f.c.water.site30)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.water.site30)$field_season),
       col=col)
anova(water.bd.site30)

phy97.f.c.water.site35 <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$site==35, prune=TRUE)
water.ufdist.site35 <- UniFrac(phy97.f.c.water.site35, weighted=T, normalized=F, parallel=F, fast=T)
water.bd.site35 <- betadisper(d=water.ufdist.site35, 
                                 group=sample_data(phy97.f.c.water.site35)$field_season,
                                 type="centroid", bias.adjust=FALSE)
plot(water.bd.site35, hull=F, label=F, 
     main="water - 35", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd.site35, sample_data(phy97.f.c.water.site35)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.water.site35)$field_season),
       col=col)
anova(water.bd.site35)

phy97.f.c.water.site8 <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$site==8, prune=TRUE)
water.ufdist.site8 <- UniFrac(phy97.f.c.water.site8, weighted=T, normalized=F, parallel=F, fast=T)
water.bd.site8 <- betadisper(d=water.ufdist.site8, 
                                group=sample_data(phy97.f.c.water.site8)$field_season,
                                type="centroid", bias.adjust=FALSE)
plot(water.bd.site8, hull=F, label=F, 
     main="water - 8", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(water.bd.site8, sample_data(phy97.f.c.water.site8)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.water.site8)$field_season),
       col=col)
anova(water.bd.site8)
dev.off()

jpeg(filename="figures/water_betadisper_boxplot_bysite.jpg", width = 13.5, height = 4, units="in",res = 300)
par(mfrow=c(1,3))
boxplot(water.bd.site30,col=col,main="Site 30")
boxplot(water.bd.site35,col=col,main="Site 35")
boxplot(water.bd.site8,col=col,main="Site 8")
dev.off()



phy97.f.c.coral <- subset_samples(phy97.f.c.coral, sample_data(phy97.f.c.coral)$site != "34")
phy97.f.c.coral <- subset_samples(phy97.f.c.coral, sample_data(phy97.f.c.coral)$field_season != "KI2015a_Pre")


phy97.f.c.coral.Mfol <- subset_samples(phy97.f.c.coral, sample_data(phy97.f.c.coral)$Coral_Species == "Montipora_foliosa")
phy97.f.c.coral.Peydo <- subset_samples(phy97.f.c.coral, sample_data(phy97.f.c.coral)$Coral_Species == "Pocillopora_eydouxi")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral, sample_data(phy97.f.c.coral)$Coral_Species == "Porites_lobata")

coral.ufdist.Mfol <- UniFrac(phy97.f.c.coral.Mfol, weighted=T, normalized=F, parallel=F, fast=T)
coral.bd.Mfol <- betadisper(d=coral.ufdist.Mfol, 
                             group=sample_data(phy97.f.c.coral.Mfol)$field_season,
                             type="centroid", bias.adjust=FALSE)
plot(coral.bd.Mfol, hull=F, label=F, 
     main="Montipora", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(coral.bd.Mfol, sample_data(phy97.f.c.coral.Mfol)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.coral.Mfol)$field_season),
       col=col)
anova(coral.bd.Mfol)

coral.ufdist.Peydo <- UniFrac(phy97.f.c.coral.Peydo, weighted=T, normalized=F, parallel=F, fast=T)
coral.bd.Peydo <- betadisper(d=coral.ufdist.Peydo, 
                            group=sample_data(phy97.f.c.coral.Peydo)$field_season,
                            type="centroid", bias.adjust=FALSE)
plot(coral.bd.Peydo, hull=F, label=F, 
     main="P. eydouxi", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(coral.bd.Peydo, sample_data(phy97.f.c.coral.Peydo)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.coral.Peydo)$field_season),
       col=col)
anova(coral.bd.Peydo)

coral.ufdist.Plob <- UniFrac(phy97.f.c.coral.Plob, weighted=T, normalized=F, parallel=F, fast=T)
coral.bd.Plob <- betadisper(d=coral.ufdist.Plob, 
                            group=sample_data(phy97.f.c.coral.Plob)$field_season,
                            type="centroid", bias.adjust=FALSE)
plot(coral.bd.Plob, hull=F, label=F, 
     main="P. lobata", col=col,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(coral.bd.Plob, sample_data(phy97.f.c.coral.Plob)$field_season,  
         draw = c("polygon"),
         col = col, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(sample_data(phy97.f.c.coral.Plob)$field_season),
       col=col)
anova(coral.bd.Plob)

boxplot(coral.bd.Mfol,col=col,main="Montipora")
boxplot(coral.bd.Peydo,col=col, main="P. eydouxi")
boxplot(coral.bd.Plob,col=col, main="P. lobata")
