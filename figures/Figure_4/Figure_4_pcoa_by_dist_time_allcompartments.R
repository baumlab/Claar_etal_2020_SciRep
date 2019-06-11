# KI Compartment Betadispersion

# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Load necessary packages
library(vegan)
library(phyloseq)

# Set seed for reproducible results
set.seed(1010)

phyASV.f.c.coral.Peyd.VH <- subset_samples(phyASV.f.c.coral.p,sample_data(phyASV.f.c.coral.p)$Dist=="VeryHigh" & sample_data(phyASV.f.c.coral.p)$Coral_Species=="Pocillopora_eydouxi")
phyASV.f.c.coral.Peyd.M <- subset_samples(phyASV.f.c.coral.p,sample_data(phyASV.f.c.coral.p)$Dist=="Medium" & sample_data(phyASV.f.c.coral.p)$Coral_Species=="Pocillopora_eydouxi")

phyASV.f.c.coral.MAeq.VH <- subset_samples(phyASV.f.c.coral.p,sample_data(phyASV.f.c.coral.p)$Dist=="VeryHigh" & sample_data(phyASV.f.c.coral.p)$Coral_Species=="Montipora_foliosa")
phyASV.f.c.coral.MAeq.M <- subset_samples(phyASV.f.c.coral.p,sample_data(phyASV.f.c.coral.p)$Dist=="Medium" & sample_data(phyASV.f.c.coral.p)$Coral_Species=="Montipora_foliosa")

phyASV.f.c.coral.Plob.VH <- subset_samples(phyASV.f.c.coral.p,sample_data(phyASV.f.c.coral.p)$Dist=="VeryHigh" & sample_data(phyASV.f.c.coral.p)$Coral_Species=="Porites_lobata")
phyASV.f.c.coral.Plob.M <- subset_samples(phyASV.f.c.coral.p,sample_data(phyASV.f.c.coral.p)$Dist=="Medium" & sample_data(phyASV.f.c.coral.p)$Coral_Species=="Porites_lobata")

phyASV.f.c.sediment.VH <- subset_samples(phyASV.f.c.sediment.p,sample_data(phyASV.f.c.sediment.p)$Dist=="VeryHigh")
phyASV.f.c.sediment.M <- subset_samples(phyASV.f.c.sediment.p,sample_data(phyASV.f.c.sediment.p)$Dist=="Medium")

phyASV.f.c.water.VH <- subset_samples(phyASV.f.c.water.p,sample_data(phyASV.f.c.water.p)$Dist=="VeryHigh")
phyASV.f.c.water.M <- subset_samples(phyASV.f.c.water.p,sample_data(phyASV.f.c.water.p)$Dist=="Medium")

# Calculate unifrac distances
sediment.ufdist <- UniFrac(phyASV.f.c.sediment, weighted=T, 
                           normalized=F, parallel=F, fast=T)
water.ufdist <- UniFrac(phyASV.f.c.water, weighted=T, 
                        normalized=F, parallel=F, fast=T)
sediment.VH.ufdist <- UniFrac(phyASV.f.c.sediment.VH, weighted=T, 
                              normalized=F, parallel=F, fast=T)
water.VH.ufdist <- UniFrac(phyASV.f.c.water.VH, weighted=T, 
                           normalized=F, parallel=F, fast=T)
sediment.M.ufdist <- UniFrac(phyASV.f.c.sediment.M, weighted=T, 
                             normalized=F, parallel=F, fast=T)
water.M.ufdist <- UniFrac(phyASV.f.c.water.M, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Peyd.VH.ufdist <- UniFrac(phyASV.f.c.coral.Peyd.VH, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Plob.VH.ufdist <- UniFrac(phyASV.f.c.coral.Plob.VH, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Maeq.VH.ufdist <- UniFrac(phyASV.f.c.coral.MAeq.VH, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Peyd.M.ufdist <- UniFrac(phyASV.f.c.coral.Peyd.M, weighted=T, 
                         normalized=F, parallel=F, fast=T)
Plob.M.ufdist <- UniFrac(phyASV.f.c.coral.Plob.M, weighted=T, 
                         normalized=F, parallel=F, fast=T)
Maeq.M.ufdist <- UniFrac(phyASV.f.c.coral.MAeq.M, weighted=T, 
                         normalized=F, parallel=F, fast=T)

# Run betadisper
sediment.bd.dist <- betadisper(d=sediment.ufdist, 
                               group=sample_data(phyASV.f.c.sediment)$Dist,
                               type="centroid", bias.adjust=FALSE)
water.bd.dist <- betadisper(d=water.ufdist, 
                            group=sample_data(phyASV.f.c.water)$Dist,
                            type="centroid", bias.adjust=FALSE)
sediment.VH.bd <- betadisper(d=sediment.VH.ufdist, 
                             group=sample_data(phyASV.f.c.sediment.VH)$field_season,
                             type="centroid", bias.adjust=FALSE)
water.VH.bd <- betadisper(d=water.VH.ufdist, 
                          group=sample_data(phyASV.f.c.water.VH)$field_season,
                          type="centroid", bias.adjust=FALSE)
sediment.M.bd <- betadisper(d=sediment.M.ufdist, 
                            group=sample_data(phyASV.f.c.sediment.M)$field_season,
                            type="centroid", bias.adjust=FALSE)
water.M.bd <- betadisper(d=water.M.ufdist, 
                         group=sample_data(phyASV.f.c.water.M)$field_season,
                         type="centroid", bias.adjust=FALSE)
Peyd.VH.bd <- betadisper(d=Peyd.VH.ufdist, 
                         group=sample_data(phyASV.f.c.coral.Peyd.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Plob.VH.bd <- betadisper(d=Plob.VH.ufdist, 
                         group=sample_data(phyASV.f.c.coral.Plob.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Maeq.VH.bd <- betadisper(d=Maeq.VH.ufdist, 
                         group=sample_data(phyASV.f.c.coral.MAeq.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Peyd.M.bd <- betadisper(d=Peyd.M.ufdist, 
                        group=sample_data(phyASV.f.c.coral.Peyd.M)$field_season,
                        type="centroid", bias.adjust=FALSE)
Plob.M.bd <- betadisper(d=Plob.M.ufdist, 
                        group=sample_data(phyASV.f.c.coral.Plob.M)$field_season,
                        type="centroid", bias.adjust=FALSE)
Maeq.M.bd <- betadisper(d=Maeq.M.ufdist, 
                        group=sample_data(phyASV.f.c.coral.MAeq.M)$field_season,
                        type="centroid", bias.adjust=FALSE)

# Stats for Table 1
anova(Peyd.VH.bd)
anova(Maeq.VH.bd)
anova(Plob.VH.bd)
anova(sediment.VH.bd)
anova(water.VH.bd)

anova(Peyd.M.bd)
anova(Maeq.M.bd)
anova(Plob.M.bd)
anova(sediment.M.bd)
anova(water.M.bd)

# Extract sample data as metadata
Peyd.VH.metadata <- as(sample_data(phyASV.f.c.coral.Peyd.VH), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.Peyd.VH <- adonis(Peyd.VH.ufdist ~ field_season, data = Peyd.VH.metadata)

#Extract sample data as metadata
Maeq.VH.metadata <- as(sample_data(phyASV.f.c.coral.MAeq.VH), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.Maeq.VH <- adonis(Maeq.VH.ufdist ~ field_season, data = Maeq.VH.metadata)

#Extract sample data as metadata
Plob.VH.metadata <- as(sample_data(phyASV.f.c.coral.Plob.VH), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.Plob.VH <- adonis(Plob.VH.ufdist ~ field_season, data = Plob.VH.metadata)

# Extract sample data as metadata
sediment.VH.metadata <- as(sample_data(phyASV.f.c.sediment.VH), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.sediment.VH <- adonis(sediment.VH.ufdist ~ field_season, data = sediment.VH.metadata)

# Extract sample data as metadata
water.VH.metadata <- as(sample_data(phyASV.f.c.water.VH), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.water.VH <- adonis(water.VH.ufdist ~ field_season, data = water.VH.metadata)

## MEDIUM
# Extract sample data as metadata
Peyd.M.metadata <- as(sample_data(phyASV.f.c.coral.Peyd.M), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.Peyd.M <- adonis(Peyd.M.ufdist ~ field_season, data = Peyd.M.metadata)

# Extract sample data as metadata
Maeq.M.metadata <- as(sample_data(phyASV.f.c.coral.MAeq.M), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.Maeq.M <- adonis(Maeq.M.ufdist ~ field_season, data = Maeq.M.metadata)

# Extract sample data as metadata
Plob.M.metadata <- as(sample_data(phyASV.f.c.coral.Plob.M), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.Plob.M <- adonis(Plob.M.ufdist ~ field_season, data = Plob.M.metadata)

# Extract sample data as metadata
sediment.M.metadata <- as(sample_data(phyASV.f.c.sediment.M), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.sediment.M <- adonis(sediment.M.ufdist ~ field_season, data = sediment.M.metadata)

# Extract sample data as metadata
water.M.metadata <- as(sample_data(phyASV.f.c.water.M), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
ado.water.M <- adonis(water.M.ufdist ~ field_season, data = water.M.metadata)

ado.Peyd.M
ado.Maeq.M
ado.Plob.M
ado.sediment.M
ado.water.M

ado.Peyd.VH
ado.Maeq.VH
ado.Plob.VH
ado.sediment.VH
ado.water.VH

# Make Figure 4 as jpeg
jpeg(filename="figures/Figure_4/Figure_4_pcoa_by_dist_time_allcompartments.jpg", 
     width = 8, height = 20, units="in",res = 300)
par(mfrow=c(5,2),mar=c(1,1,5,1))

# Pocillopora
plot(Peyd.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Pocillopora grandis")), line = 1, cex.main=2.2)
ordihull(Peyd.VH.bd, sample_data(phyASV.f.c.coral.Peyd.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=c(16,17,18), col=timecols, cex=2.7, pt.cex=4,
       legend=c("Before","Post Storm","After"))
mtext("VERY HIGH DISTURBANCE", side=3, line=3, cex=1.5)

plot(Peyd.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Pocillopora grandis")), line = 1,cex.main=2.2)
ordihull(Peyd.M.bd, sample_data(phyASV.f.c.coral.Peyd.M)$field_season,
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
mtext("MEDIUM DISTURBANCE", side=3, line=3, cex=1.5)

par(mar=c(1,1,2.3,1))
plot(Plob.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Porites lobata")), line = 1,cex.main=2.2)
ordihull(Plob.VH.bd, sample_data(phyASV.f.c.coral.Plob.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(Plob.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Porites lobata")), line = 1, cex.main=2.2)
ordihull(Plob.M.bd, sample_data(phyASV.f.c.coral.Plob.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(Maeq.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Montipora aequituberculata")), line = 1, cex.main=2.2)
ordihull(Maeq.VH.bd, sample_data(phyASV.f.c.coral.MAeq.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(Maeq.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Montipora aequituberculata")), line = 1, cex.main=2.2)
ordihull(Maeq.M.bd, sample_data(phyASV.f.c.coral.MAeq.M)$field_season, 
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

# Sediment
plot(sediment.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Sediment", line = 1, cex.main=2.2, font.main=1)
ordihull(sediment.VH.bd, sample_data(phyASV.f.c.sediment.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(sediment.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Sediment", line = 1, cex.main=2.2, font.main=1)
ordihull(sediment.M.bd, sample_data(phyASV.f.c.sediment.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

# Water
plot(water.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Water", line = 1, cex.main=2.2, font.main=1)
ordihull(water.VH.bd, sample_data(phyASV.f.c.water.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(water.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Water", line = 1, cex.main=2.2, font.main=1)
ordihull(water.M.bd, sample_data(phyASV.f.c.water.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

dev.off()

# Make Figure 4 as pdf
pdf(file="figures/Figure_4/Figure_4_pcoa_by_dist_time_allcompartments.pdf", 
     width = 8, height = 20, useDingbats = FALSE)
par(mfrow=c(5,2),mar=c(1,1,5,1))

# Pocillopora
plot(Peyd.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Pocillopora grandis")), line = 1, cex.main=2.2)
ordihull(Peyd.VH.bd, sample_data(phyASV.f.c.coral.Peyd.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=c(16,17,18), col=timecols, cex=2.7, pt.cex=4,
       legend=c("Before","Post Storm","After"))
mtext("VERY HIGH DISTURBANCE", side=3, line=3, cex=1.5)

plot(Peyd.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Pocillopora grandis")), line = 1,cex.main=2.2)
ordihull(Peyd.M.bd, sample_data(phyASV.f.c.coral.Peyd.M)$field_season,
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
mtext("MEDIUM DISTURBANCE", side=3, line=3, cex=1.5)

par(mar=c(1,1,2.3,1))
plot(Plob.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Porites lobata")), line = 1,cex.main=2.2)
ordihull(Plob.VH.bd, sample_data(phyASV.f.c.coral.Plob.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(Plob.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Porites lobata")), line = 1, cex.main=2.2)
ordihull(Plob.M.bd, sample_data(phyASV.f.c.coral.Plob.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(Maeq.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Montipora aequituberculata")), line = 1, cex.main=2.2)
ordihull(Maeq.VH.bd, sample_data(phyASV.f.c.coral.MAeq.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(Maeq.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Montipora aequituberculata")), line = 1, cex.main=2.2)
ordihull(Maeq.M.bd, sample_data(phyASV.f.c.coral.MAeq.M)$field_season, 
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

# Sediment
plot(sediment.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Sediment", line = 1, cex.main=2.2, font.main=1)
ordihull(sediment.VH.bd, sample_data(phyASV.f.c.sediment.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(sediment.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Sediment", line = 1, cex.main=2.2, font.main=1)
ordihull(sediment.M.bd, sample_data(phyASV.f.c.sediment.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

# Water
plot(water.VH.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Water", line = 1, cex.main=2.2, font.main=1)
ordihull(water.VH.bd, sample_data(phyASV.f.c.water.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

plot(water.M.bd, hull=F, label=F, cex=2,
     main=NULL, col=timecols, pch=c(16,17,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Water", line = 1, cex.main=2.2, font.main=1)
ordihull(water.M.bd, sample_data(phyASV.f.c.water.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)

dev.off()

