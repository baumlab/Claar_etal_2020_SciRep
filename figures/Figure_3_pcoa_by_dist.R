# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)
library(ggplot2)

phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")


sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, 
                       normalized=F, parallel=F, fast=T)
sediment.bd.dist <- betadisper(d=sediment.ufdist, 
                           group=sample_data(phy97.f.c.sediment)$Dist,
                           type="centroid", bias.adjust=FALSE)

water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, 
                           normalized=F, parallel=F, fast=T)
water.bd.dist <- betadisper(d=water.ufdist, 
                               group=sample_data(phy97.f.c.water)$Dist,
                               type="centroid", bias.adjust=FALSE)

all.ufdist <- UniFrac(phy97.f.c, weighted=T, 
                        normalized=F, parallel=F, fast=T)
all.bd.dist <- betadisper(d=all.ufdist, 
                            group=sample_data(phy97.f.c)$Dist,
                            type="centroid", bias.adjust=FALSE)


Peyd.ufdist <- UniFrac(phy97.f.c.coral.Peyd, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Peyd.bd.dist <- betadisper(d=Peyd.ufdist, 
                           group=sample_data(phy97.f.c.coral.Peyd)$Dist,
                           type="centroid", bias.adjust=FALSE)

MAeq.ufdist <- UniFrac(phy97.f.c.coral.MAeq, weighted=T, 
                       normalized=F, parallel=F, fast=T)
MAeq.bd.dist <- betadisper(d=MAeq.ufdist, 
                           group=sample_data(phy97.f.c.coral.MAeq)$Dist,
                           type="centroid", bias.adjust=FALSE)

Plob.ufdist <- UniFrac(phy97.f.c.coral.Plob, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Plob.bd.dist <- betadisper(d=Plob.ufdist, 
                           group=sample_data(phy97.f.c.coral.Plob)$Dist,
                           type="centroid", bias.adjust=FALSE)

coral.ufdist <- UniFrac(phy97.f.c.coral, weighted=T, 
                        normalized=F, parallel=F, fast=T)
coral.bd.dist <- betadisper(d=coral.ufdist, 
                            group=sample_data(phy97.f.c.coral)$Dist,
                            type="centroid", bias.adjust=FALSE)
coral.bd.species <- betadisper(d=coral.ufdist, 
                               group=sample_data(phy97.f.c.coral)$Coral_Species,
                               type="centroid", bias.adjust=FALSE)

betadisper.sediment <- anova(sediment.bd.dist)
betadisper.water <- anova(water.bd.dist)
betadisper.Peyd <- anova(Peyd.bd.dist)
betadisper.MAeq <- anova(MAeq.bd.dist)
betadisper.Plob <- anova(Plob.bd.dist)
betadisper.all <- anova(all.bd.dist)

betadisper.sediment$`Pr(>F)`[1]
betadisper.sediment$`F value`[1]

set.seed(2020)
sediment.metadata <- as(sample_data(phy97.f.c.sediment), "data.frame")
adonis(distance(phy97.f.c.sediment, method="wunifrac") ~ Dist,
       data = sediment.metadata)

water.metadata <- as(sample_data(phy97.f.c.water), "data.frame")
adonis(distance(phy97.f.c.water, method="wunifrac") ~ Dist,
       data = water.metadata)

plob.metadata <- as(sample_data(phy97.f.c.coral.Plob), "data.frame")
adonis(distance(phy97.f.c.coral.Plob, method="wunifrac") ~ Dist,
       data = plob.metadata)

peyd.metadata <- as(sample_data(phy97.f.c.coral.Peyd), "data.frame")
adonis(distance(phy97.f.c.coral.Peyd, method="wunifrac") ~ Dist,
       data = peyd.metadata)

maeq.metadata <- as(sample_data(phy97.f.c.coral.MAeq), "data.frame")
adonis(distance(phy97.f.c.coral.MAeq, method="wunifrac") ~ Dist,
       data = maeq.metadata)

all.metadata <- as(sample_data(phy97.f.c), "data.frame")
adonis(distance(phy97.f.c, method="wunifrac") ~ Dist,
       data = all.metadata)


# jpeg(filename="figures/pcoa_by_dist2.jpg", 
jpeg(filename="figures/Figure_3_pcoa_by_dist.jpg", 
     width = 12, height = 8, units="in",res = 300)
par(mfrow=c(2,3), mar=c(1,1,3,1))
# Pocillopora
plot(Peyd.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Pocillopora grandis"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Peyd.bd.dist, sample_data(phy97.f.c.coral.Peyd)$Dist,  
         draw = c("polygon"), col = sitecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=1:6, cex=2.5,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd)$Dist)),
       col=sitecols)
text(-0.07,-0.037,"betadisper = sig. ***",cex=2)
text(-0.07,-0.046,"adonis = sig. ***",cex=2)
mtext("A)", side=3,line=0.5,adj=-0.025,cex=2) 


# Montipora
plot(MAeq.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("  Montipora aequituburculata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(MAeq.bd.dist, sample_data(phy97.f.c.coral.MAeq)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq)$Dist)),
#        col=sitecols)
text(-0.0006,-0.00188,"betadisper = sig. ***",cex=2)
text(-0.0006,-0.00235,"adonis = sig. *",cex=2)
mtext("B)", side=3,line=0.5,adj=-0.025,cex=2) 


# Porites
plot(Plob.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Porites lobata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Plob.bd.dist, sample_data(phy97.f.c.coral.Plob)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob)$Dist)),
#        col=sitecols)
text(-0.00275,-0.018,"betadisper = not sig.",cex=2)
text(-0.00275,-0.02,"adonis = not sig.",cex=2)
mtext("C)", side=3,line=0.5,adj=-0.025,cex=2) 


# Sediment
plot(sediment.bd.dist, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(sediment.bd.dist, sample_data(phy97.f.c.sediment)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.sediment)$Dist)),
#        col=sitecols)
text(-0.02,-0.075,"betadisper = sig. ***",cex=2)
text(-0.02,-0.087,"adonis = sig. ***",cex=2)
mtext("D)", side=3,line=0.5,adj=-0.025,cex=2) 


# Water
plot(water.bd.dist, hull=F, label=F, 
     main="Water", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(water.bd.dist, sample_data(phy97.f.c.water)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$Dist)),
#        col=sitecols)
text(0.045,-0.040,"betadisper = not sig.",cex=2)
text(0.045,-0.047,"adonis = sig. **",cex=2)
mtext("E)", side=3,line=0.5,adj=-0.025,cex=2) 


# All
plot(all.bd.dist, hull=F, label=F, 
     main="All Compartments", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(all.bd.dist, sample_data(phy97.f.c)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$Dist)),
#        col=sitecols)
text(0.065,-0.027,"betadisper = sig. ***",cex=2)
text(0.065,-0.037,"adonis = sig. ***",cex=2) 
mtext("F)", side=3,line=0.5,adj=-0.025,cex=2) 


dev.off()



pdf(file="figures/Figure_3_pcoa_by_dist.pdf", 
     width = 12, height = 8,useDingbats = FALSE)
par(mfrow=c(2,3), mar=c(1,1,3,1))
# Pocillopora
plot(Peyd.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Pocillopora grandis"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Peyd.bd.dist, sample_data(phy97.f.c.coral.Peyd)$Dist,  
         draw = c("polygon"), col = sitecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=1:6, cex=2.5,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd)$Dist)),
       col=sitecols)
text(-0.07,-0.037,"betadisper = sig. ***",cex=2)
text(-0.07,-0.046,"adonis = sig. ***",cex=2)
mtext("A)", side=3,line=0.5,adj=-0.025,cex=2) 


# Montipora
plot(MAeq.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("  Montipora aequituburculata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(MAeq.bd.dist, sample_data(phy97.f.c.coral.MAeq)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq)$Dist)),
#        col=sitecols)
text(-0.0006,-0.00188,"betadisper = sig. ***",cex=2)
text(-0.0006,-0.00235,"adonis = sig. *",cex=2)
mtext("B)", side=3,line=0.5,adj=-0.025,cex=2) 


# Porites
plot(Plob.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Porites lobata"))), col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(Plob.bd.dist, sample_data(phy97.f.c.coral.Plob)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob)$Dist)),
#        col=sitecols)
text(-0.00275,-0.018,"betadisper = not sig.",cex=2)
text(-0.00275,-0.02,"adonis = not sig.",cex=2)
mtext("C)", side=3,line=0.5,adj=-0.025,cex=2) 


# Sediment
plot(sediment.bd.dist, hull=F, label=F, 
     main="Sediment", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(sediment.bd.dist, sample_data(phy97.f.c.sediment)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.sediment)$Dist)),
#        col=sitecols)
text(-0.02,-0.075,"betadisper = sig. ***",cex=2)
text(-0.02,-0.087,"adonis = sig. ***",cex=2)
mtext("D)", side=3,line=0.5,adj=-0.025,cex=2) 


# Water
plot(water.bd.dist, hull=F, label=F, 
     main="Water", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(water.bd.dist, sample_data(phy97.f.c.water)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$Dist)),
#        col=sitecols)
text(0.045,-0.040,"betadisper = not sig.",cex=2)
text(0.045,-0.047,"adonis = sig. **",cex=2)
mtext("E)", side=3,line=0.5,adj=-0.025,cex=2) 


# All
plot(all.bd.dist, hull=F, label=F, 
     main="All Compartments", col=sitecols,
     xlab="", ylab="", cex=2, sub="", 
     xaxt='n',yaxt='n',cex.main=2.5)
ordihull(all.bd.dist, sample_data(phy97.f.c)$Dist, 
         draw = c("polygon"),
         col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6, cex=2.5,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$Dist)),
#        col=sitecols)
text(0.065,-0.027,"betadisper = sig. ***",cex=2)
text(0.065,-0.037,"adonis = sig. ***",cex=2) 
mtext("F)", side=3,line=0.5,adj=-0.025,cex=2) 


dev.off()