# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)
library(ggplot2)
library(phyloseq)

phyASV.f.c.coral.Peyd <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phyASV.f.c.coral.MAeq <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Montipora_foliosa")
phyASV.f.c.coral.Plob <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Porites_lobata")

sediment.ufdist <- UniFrac(phyASV.f.c.sediment, weighted=T, 
                       normalized=F, parallel=F, fast=T)
sediment.bd.dist <- betadisper(d=sediment.ufdist, 
                           group=sample_data(phyASV.f.c.sediment)$Dist,
                           type="centroid", bias.adjust=FALSE)

water.ufdist <- UniFrac(phyASV.f.c.water, weighted=T, 
                           normalized=F, parallel=F, fast=T)
water.bd.dist <- betadisper(d=water.ufdist, 
                               group=sample_data(phyASV.f.c.water)$Dist,
                               type="centroid", bias.adjust=FALSE)

all.ufdist <- UniFrac(phyASV.f.c, weighted=T, 
                        normalized=F, parallel=F, fast=T)
all.bd.dist <- betadisper(d=all.ufdist, 
                            group=sample_data(phyASV.f.c)$Dist,
                            type="centroid", bias.adjust=FALSE)

st.uf <- UniFrac(phyASV.f.c, weighted=T, 
                      normalized=F, parallel=F, fast=T)
st.bd <- betadisper(d=st.uf, 
                          group=sample_data(phyASV.f.c)$SampleType,
                          type="centroid", bias.adjust=FALSE)

Peyd.ufdist <- UniFrac(phyASV.f.c.coral.Peyd, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Peyd.bd.dist <- betadisper(d=Peyd.ufdist, 
                           group=sample_data(phyASV.f.c.coral.Peyd)$Dist,
                           type="centroid", bias.adjust=FALSE)

MAeq.ufdist <- UniFrac(phyASV.f.c.coral.MAeq, weighted=T, 
                       normalized=F, parallel=F, fast=T)
MAeq.bd.dist <- betadisper(d=MAeq.ufdist, 
                           group=sample_data(phyASV.f.c.coral.MAeq)$Dist,
                           type="centroid", bias.adjust=FALSE)

Plob.ufdist <- UniFrac(phyASV.f.c.coral.Plob, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Plob.bd.dist <- betadisper(d=Plob.ufdist, 
                           group=sample_data(phyASV.f.c.coral.Plob)$Dist,
                           type="centroid", bias.adjust=FALSE)

coral.ufdist <- UniFrac(phyASV.f.c.coral, weighted=T, 
                        normalized=F, parallel=F, fast=T)
coral.bd.dist <- betadisper(d=coral.ufdist, 
                            group=sample_data(phyASV.f.c.coral)$Dist,
                            type="centroid", bias.adjust=FALSE)
coral.bd.species <- betadisper(d=coral.ufdist, 
                               group=sample_data(phyASV.f.c.coral)$Coral_Species,
                               type="centroid", bias.adjust=FALSE)

betadisper.sediment <- anova(sediment.bd.dist)
betadisper.water <- anova(water.bd.dist)
betadisper.Peyd <- anova(Peyd.bd.dist)
betadisper.MAeq <- anova(MAeq.bd.dist)
betadisper.Plob <- anova(Plob.bd.dist)
betadisper.all <- anova(all.bd.dist)

betadisper.sediment
betadisper.water
betadisper.Peyd
betadisper.MAeq
betadisper.Plob
betadisper.all

betadisper.sediment$`Pr(>F)`[1]
betadisper.sediment$`F value`[1]

set.seed(2020)
sediment.metadata <- as(sample_data(phyASV.f.c.sediment), "data.frame")
ado.sediment <- adonis(sediment.ufdist ~ Dist, data = sediment.metadata)

water.metadata <- as(sample_data(phyASV.f.c.water), "data.frame")
ado.water <- adonis(water.ufdist ~ Dist, data = water.metadata)

Plob.metadata <- as(sample_data(phyASV.f.c.coral.Plob), "data.frame")
ado.Plob <- adonis(Plob.ufdist ~ Dist, data = Plob.metadata)

Peyd.metadata <- as(sample_data(phyASV.f.c.coral.Peyd), "data.frame")
ado.Peyd <- adonis(Peyd.ufdist ~ Dist, data = Peyd.metadata)

Maeq.metadata <- as(sample_data(phyASV.f.c.coral.MAeq), "data.frame")
ado.Maeq <- adonis(MAeq.ufdist ~ Dist, data = Maeq.metadata)

all.metadata <- as(sample_data(phyASV.f.c), "data.frame")
ado.all <- adonis(all.ufdist ~ Dist, data = all.metadata)

ado.sediment
ado.water
ado.Peyd
ado.Maeq
ado.Plob
ado.all

# Make figure 3 as jpeg
jpeg(filename="figures/Figure_3/Figure_3_pcoa_by_dist.jpg", 
     width = 12, height = 8, units="in",res = 300)
par(mfrow=c(2,3), mar=c(1,1,3,1))
# Pocillopora
plot(Peyd.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Pocillopora grandis"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Peyd.bd.dist, sample_data(phyASV.f.c.coral.Peyd)$Dist,  
         draw = c("polygon"), col = sitecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=1:6, cex=2.5,
       legend=levels(as.factor(sample_data(phyASV.f.c.coral.Peyd)$Dist)),
       col=sitecols)
text(0.023,-0.044,"betadisper = sig. ***",cex=2)
text(0.023,-0.052,"adonis = sig. ***",cex=2)
mtext("A)", side=3,line=0.5,adj=-0.025,cex=2) 


# Montipora
plot(MAeq.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("  Montipora aequituburculata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(MAeq.bd.dist, sample_data(phyASV.f.c.coral.MAeq)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.0012,-0.0011,"betadisper = sig. **",cex=2)
text(0.0012,-0.0015,"adonis = not sig.",cex=2)
mtext("B)", side=3,line=0.5,adj=-0.025,cex=2) 


# Porites
plot(Plob.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Porites lobata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Plob.bd.dist, sample_data(phyASV.f.c.coral.Plob)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(-0.00415,0.024,"betadisper = not sig.",cex=2)
text(-0.00415,0.0225,"adonis = not sig.",cex=2)
mtext("C)", side=3,line=0.5,adj=-0.025,cex=2) 


# Sediment
plot(sediment.bd.dist, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(sediment.bd.dist, sample_data(phyASV.f.c.sediment)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.076,0.112,"betadisper = sig. ***",cex=2)
text(0.076,0.1,"adonis = sig. ***",cex=2)
mtext("D)", side=3,line=0.5,adj=-0.025,cex=2) 


# Water
plot(water.bd.dist, hull=F, label=F, 
     main="Water", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(water.bd.dist, sample_data(phyASV.f.c.water)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.08,0.09,"betadisper = sig. **",cex=2)
text(0.08,0.08,"adonis = sig. **",cex=2)
mtext("E)", side=3,line=0.5,adj=-0.025,cex=2) 


# All
plot(all.bd.dist, hull=F, label=F, 
     main="All Compartments", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(all.bd.dist, sample_data(phyASV.f.c)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.017,0.089,"betadisper = sig. ***",cex=2)
text(0.017,0.077,"adonis = sig. ***",cex=2) 
mtext("F)", side=3,line=0.5,adj=-0.025,cex=2) 

dev.off()

# Make same figure as pdf
pdf(file="figures/Figure_3/Figure_3_pcoa_by_dist.pdf",
     width = 12, height = 8,useDingbats = FALSE)
par(mfrow=c(2,3), mar=c(1,1,3,1))
# Pocillopora
plot(Peyd.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Pocillopora grandis"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Peyd.bd.dist, sample_data(phyASV.f.c.coral.Peyd)$Dist,  
         draw = c("polygon"), col = sitecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=1:6, cex=2.5,
       legend=levels(as.factor(sample_data(phyASV.f.c.coral.Peyd)$Dist)),
       col=sitecols)
text(0.023,-0.044,"betadisper = sig. ***",cex=2)
text(0.023,-0.052,"adonis = sig. ***",cex=2)
mtext("A)", side=3,line=0.5,adj=-0.025,cex=2) 


# Montipora
plot(MAeq.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("  Montipora aequituburculata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(MAeq.bd.dist, sample_data(phyASV.f.c.coral.MAeq)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.0012,-0.0011,"betadisper = sig. **",cex=2)
text(0.0012,-0.0015,"adonis = not sig.",cex=2)
mtext("B)", side=3,line=0.5,adj=-0.025,cex=2) 


# Porites
plot(Plob.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Porites lobata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Plob.bd.dist, sample_data(phyASV.f.c.coral.Plob)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(-0.00415,0.024,"betadisper = not sig.",cex=2)
text(-0.00415,0.0225,"adonis = not sig.",cex=2)
mtext("C)", side=3,line=0.5,adj=-0.025,cex=2) 


# Sediment
plot(sediment.bd.dist, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(sediment.bd.dist, sample_data(phyASV.f.c.sediment)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.076,0.112,"betadisper = sig. ***",cex=2)
text(0.076,0.1,"adonis = sig. ***",cex=2)
mtext("D)", side=3,line=0.5,adj=-0.025,cex=2) 


# Water
plot(water.bd.dist, hull=F, label=F, 
     main="Water", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(water.bd.dist, sample_data(phyASV.f.c.water)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.08,0.09,"betadisper = sig. **",cex=2)
text(0.08,0.08,"adonis = sig. **",cex=2)
mtext("E)", side=3,line=0.5,adj=-0.025,cex=2) 


# All
plot(all.bd.dist, hull=F, label=F, 
     main="All Compartments", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(all.bd.dist, sample_data(phyASV.f.c)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
text(0.017,0.089,"betadisper = sig. ***",cex=2)
text(0.017,0.077,"adonis = sig. ***",cex=2) 
mtext("F)", side=3,line=0.5,adj=-0.025,cex=2) 
dev.off()

