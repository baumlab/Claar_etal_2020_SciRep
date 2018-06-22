# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)

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


jpeg(filename="figures/pcoa_by_dist_time_coralsp.jpg", 
     width = 8, height = 12, units="in",res = 300)
par(mfrow=c(3,2),mar=c(4,5,5,2))



plot(Peyd.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
                xlab="PCoA 1", ylab="PCoA 2", sub="")
title(expression(italic("Pocillopora grandis")), line = 1)
ordihull(Peyd.VH.bd, sample_data(phy97.f.c.coral.Peyd.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd.VH)$field_season)),
       col=timecols)
mtext("VERY HIGH DISTURBANCE",side=3,line=3)

plot(Peyd.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title(expression(italic("Pocillopora grandis")), line = 1)
ordihull(Peyd.M.bd, sample_data(phy97.f.c.coral.Peyd.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd.VH)$field_season)),
       col=timecols)
mtext("MEDIUM DISTURBANCE",side=3,line=3)

plot(Plob.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title(expression(italic("Porites lobata")), line = 1)
ordihull(Plob.VH.bd, sample_data(phy97.f.c.coral.Plob.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob.VH)$field_season)),
       col=timecols)

plot(Plob.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="",select=c(""))
title(expression(italic("Porites lobata")), line = 1)
ordihull(Plob.M.bd, sample_data(phy97.f.c.coral.Plob.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob.VH)$field_season)),
       col=timecols)

plot(Maeq.VH.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title(expression(italic("Montipora aequituberculata")), line = 1)
ordihull(Maeq.VH.bd, sample_data(phy97.f.c.coral.MAeq.VH)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq.VH)$field_season)),
       col=timecols)

plot(Maeq.M.bd, hull=F, label=F, 
     main=NULL, col=timecols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
title(expression(italic("Montipora aequituberculata")), line = 1)
ordihull(Maeq.M.bd, sample_data(phy97.f.c.coral.MAeq.M)$field_season,  draw = c("polygon"),
         col = timecols, alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=1:6,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq.VH)$field_season)),
       col=timecols)

dev.off()