# KI Compartment Betadispersion

# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Load necessary packages
library(vegan)
library(phyloseq)

# Calculate unifrac distances
sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, 
                           normalized=F, parallel=F, fast=T)
water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, 
                        normalized=F, parallel=F, fast=T)
sediment.VH.ufdist <- UniFrac(phy97.f.c.sediment.VH, weighted=T, 
                              normalized=F, parallel=F, fast=T)
water.VH.ufdist <- UniFrac(phy97.f.c.water.VH, weighted=T, 
                           normalized=F, parallel=F, fast=T)
sediment.M.ufdist <- UniFrac(phy97.f.c.sediment.M, weighted=T, 
                             normalized=F, parallel=F, fast=T)
water.M.ufdist <- UniFrac(phy97.f.c.water.M, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Peyd.VH.ufdist <- UniFrac(phy97.f.c.coral.Peyd.VH, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Plob.VH.ufdist <- UniFrac(phy97.f.c.coral.Plob.VH, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Maeq.VH.ufdist <- UniFrac(phy97.f.c.coral.MAeq.VH, weighted=T, 
                          normalized=F, parallel=F, fast=T)
Peyd.M.ufdist <- UniFrac(phy97.f.c.coral.Peyd.M, weighted=T, 
                         normalized=F, parallel=F, fast=T)
Plob.M.ufdist <- UniFrac(phy97.f.c.coral.Plob.M, weighted=T, 
                         normalized=F, parallel=F, fast=T)
Maeq.M.ufdist <- UniFrac(phy97.f.c.coral.MAeq.M, weighted=T, 
                         normalized=F, parallel=F, fast=T)

# Run betadisper
sediment.bd.dist <- betadisper(d=sediment.ufdist, 
                               group=sample_data(phy97.f.c.sediment)$Dist,
                               type="centroid", bias.adjust=FALSE)
water.bd.dist <- betadisper(d=water.ufdist, 
                            group=sample_data(phy97.f.c.water)$Dist,
                            type="centroid", bias.adjust=FALSE)
sediment.VH.bd <- betadisper(d=sediment.VH.ufdist, 
                             group=sample_data(phy97.f.c.sediment.VH)$field_season,
                             type="centroid", bias.adjust=FALSE)
water.VH.bd <- betadisper(d=water.VH.ufdist, 
                          group=sample_data(phy97.f.c.water.VH)$field_season,
                          type="centroid", bias.adjust=FALSE)
sediment.M.bd <- betadisper(d=sediment.M.ufdist, 
                            group=sample_data(phy97.f.c.sediment.M)$field_season,
                            type="centroid", bias.adjust=FALSE)
water.M.bd <- betadisper(d=water.M.ufdist, 
                         group=sample_data(phy97.f.c.water.M)$field_season,
                         type="centroid", bias.adjust=FALSE)
Peyd.VH.bd <- betadisper(d=Peyd.VH.ufdist, 
                         group=sample_data(phy97.f.c.coral.Peyd.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Plob.VH.bd <- betadisper(d=Plob.VH.ufdist, 
                         group=sample_data(phy97.f.c.coral.Plob.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Maeq.VH.bd <- betadisper(d=Maeq.VH.ufdist, 
                         group=sample_data(phy97.f.c.coral.MAeq.VH)$field_season,
                         type="centroid", bias.adjust=FALSE)
Peyd.M.bd <- betadisper(d=Peyd.M.ufdist, 
                        group=sample_data(phy97.f.c.coral.Peyd.M)$field_season,
                        type="centroid", bias.adjust=FALSE)
Plob.M.bd <- betadisper(d=Plob.M.ufdist, 
                        group=sample_data(phy97.f.c.coral.Plob.M)$field_season,
                        type="centroid", bias.adjust=FALSE)
Maeq.M.bd <- betadisper(d=Maeq.M.ufdist, 
                        group=sample_data(phy97.f.c.coral.MAeq.M)$field_season,
                        type="centroid", bias.adjust=FALSE)




jpeg(filename="figures/Figure_4_pcoa_by_dist_time_allcompartments.jpg", 
     width = 8, height = 20, units="in",res = 300)
par(mfrow=c(5,2),mar=c(1,1,5,1))

# Pocillopora
plot(Peyd.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Pocillopora grandis")), line = 1,cex.main=2.2)
ordihull(Peyd.VH.bd, sample_data(phy97.f.c.coral.Peyd.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
#legend("topleft", bty="n", pch=1:6,col=timecols,
#       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd.VH)$field_season)))
mtext("VERY HIGH DISTURBANCE",side=3,line=3,cex=1.5)

plot(Peyd.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Pocillopora grandis")), line = 1,cex.main=2.2)
ordihull(Peyd.M.bd, sample_data(phy97.f.c.coral.Peyd.M)$field_season,
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=1:6, col=timecols, cex=2.7,
       legend=c("Before","Post Storm", "After"))
       # legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd.VH)$field_season)))
mtext("MEDIUM DISTURBANCE",side=3,line=3,cex=1.5)

par(mar=c(1,1,2.3,1))
plot(Plob.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Porites lobata")), line = 1,cex.main=2.2)
ordihull(Plob.VH.bd, sample_data(phy97.f.c.coral.Plob.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob.VH)$field_season)))

plot(Plob.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Porites lobata")), line = 1,cex.main=2.2)
ordihull(Plob.M.bd, sample_data(phy97.f.c.coral.Plob.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob.VH)$field_season)))

plot(Maeq.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Montipora aequituberculata")), line = 1,cex.main=2.2)
ordihull(Maeq.VH.bd, sample_data(phy97.f.c.coral.MAeq.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq.VH)$field_season)))

plot(Maeq.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("Montipora aequituberculata")), line = 1,cex.main=2.2)
ordihull(Maeq.M.bd, sample_data(phy97.f.c.coral.MAeq.M)$field_season, 
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq.VH)$field_season)))

# Sediment
plot(sediment.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols, 
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Sediment", line = 1,cex.main=2.2,font.main=1)
ordihull(sediment.VH.bd, sample_data(phy97.f.c.sediment.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("bottomright", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.sediment)$field_season)))

plot(sediment.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Sediment", line = 1,cex.main=2.2,font.main=1)
ordihull(sediment.M.bd, sample_data(phy97.f.c.sediment.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.sediment)$field_season)))

# Water
plot(water.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Water", line = 1,cex.main=2.2,font.main=1)
ordihull(water.VH.bd, sample_data(phy97.f.c.water.VH)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("bottomright", bty="n", pch=1:6, col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$field_season)))

plot(water.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("Water", line = 1,cex.main=2.2,font.main=1)
ordihull(water.M.bd, sample_data(phy97.f.c.water.M)$field_season,  
         draw = c("polygon"), col = timecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,col=timecols,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$field_season)))

dev.off()


# pdf(file="figures/Fig_5.pdf", 
#     width = 8, height = 20,useDingbats = FALSE)
# par(mfrow=c(2,5))
# 
# 
# dev.off()