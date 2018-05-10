# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)


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

jpeg(filename="figures/pcoa_by_dist.jpg", 
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
text(-0.089,-0.037,paste("F = ", round((betadisper.Peyd$`F value`[1]),
                                       digits=1)),cex=2)
text(-0.089,-0.046,paste("p ",(ifelse(betadisper.Peyd$`Pr(>F)`[1]<0.001,"< 0.001", paste("= ",round(betadisper.Peyd$`Pr(>F)`[1],digits=3))))),cex=2)
mtext("a)", side=3,line=0.5,adj=-0.025,cex=2) 


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
text(0.0001,-0.00188,paste("F = ", round((betadisper.MAeq$`F value`[1]),
                                         digits=1)),cex=2)
text(0.0001,-0.00235,paste("p ",(ifelse(betadisper.MAeq$`Pr(>F)`[1]<0.001,"< 0.001", paste("= ",round(betadisper.MAeq$`Pr(>F)`[1],digits=3))))),cex=2)
mtext("b)", side=3,line=0.5,adj=-0.025,cex=2) 


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
text(0.015,-0.02,"p > 0.05",cex=2)
mtext("c)", side=3,line=0.5,adj=-0.025,cex=2) 


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
text(-0.037,-0.075,paste("F = ", round((betadisper.sediment$`F value`[1]),
                                       digits=1)),cex=2)
text(-0.037,-0.087,paste("p ",(ifelse(betadisper.sediment$`Pr(>F)`[1]<0.001,"< 0.001", paste("= ",round(betadisper.sediment$`Pr(>F)`[1],digits=3))))),cex=2)
mtext("d)", side=3,line=0.5,adj=-0.025,cex=2) 


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
text(0.06,-0.047,"p > 0.05",cex=2)
mtext("e)", side=3,line=0.5,adj=-0.025,cex=2) 


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
text(0.086,-0.027,paste("F = ", round((betadisper.all$`F value`[1]),
                                      digits=1)),cex=2)
text(0.086,-0.037,paste("p ",(ifelse(betadisper.all$`Pr(>F)`[1]<0.001,"< 0.001", paste("= ",round(betadisper.all$`Pr(>F)`[1],digits=3))))),cex=2)
mtext("f)", side=3,line=0.5,adj=-0.025,cex=2) 


dev.off()