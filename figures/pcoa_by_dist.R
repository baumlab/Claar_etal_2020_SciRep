# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(vegan)

sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, 
                           normalized=F, parallel=F, fast=T)
water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, 
                        normalized=F, parallel=F, fast=T)
sediment.bd.dist <- betadisper(d=sediment.ufdist, 
                               group=sample_data(phy97.f.c.sediment)$Dist,
                               type="centroid", bias.adjust=FALSE)
water.bd.dist <- betadisper(d=water.ufdist, 
                            group=sample_data(phy97.f.c.water)$Dist,
                            type="centroid", bias.adjust=FALSE)

all.ufdist <- UniFrac(phy97.f.c, weighted=T, 
                      normalized=F, parallel=F, fast=T)
all.bd.dist <- betadisper(d=all.ufdist, 
                            group=sample_data(phy97.f.c)$Dist,
                            type="centroid", bias.adjust=FALSE)


phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

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


# jpeg(filename="figures/pcoa_by_dist2.jpg", 
#      width = 20, height = 8, units="in",res = 300)
# par(mfrow=c(2,5))
# 
# # Sediment
# plot(sediment.bd.dist, hull=F, label=F, 
#      main="Sediment", col=sitecols,
#      xlab="PCoA 1", ylab="PCoA 2", sub="")
# ordihull(sediment.bd.dist, sample_data(phy97.f.c.sediment)$Dist,  draw = c("polygon"),
#          col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(as.factor(sample_data(phy97.f.c.sediment)$Dist)),
#        col=sitecols)
# 
# 
# # Water
# plot(water.bd.dist, hull=F, label=F, 
#      main="Water", col=sitecols,
#      xlab="PCoA 1", ylab="PCoA 2", sub="")
# ordihull(water.bd.dist, sample_data(phy97.f.c.water)$Dist,  draw = c("polygon"),
#          col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(as.factor(sample_data(phy97.f.c.water)$Dist)),
#        col=sitecols)
# 
# # Pocillopora
# plot(Peyd.bd.dist, hull=F, label=F, 
#      main="Pocillopora eydouxi", col=sitecols,
#      xlab="PCoA 1", ylab="PCoA 2", sub="")
# ordihull(Peyd.bd.dist, sample_data(phy97.f.c.coral.Peyd)$Dist,  draw = c("polygon"),
#          col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd)$Dist)),
#        col=sitecols)
# 
# 
# # Montipora
# plot(MAeq.bd.dist, hull=F, label=F, 
#      main="Montipora aequituburculata", col=sitecols,
#      xlab="PCoA 1", ylab="PCoA 2", sub="")
# ordihull(MAeq.bd.dist, sample_data(phy97.f.c.coral.MAeq)$Dist,  draw = c("polygon"),
#          col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.MAeq)$Dist)),
#        col=sitecols)
# 
# # Porites
# plot(Plob.bd.dist, hull=F, label=F, 
#      main="Porites lobata", col=sitecols,
#      xlab="PCoA 1", ylab="PCoA 2", sub="")
# ordihull(Plob.bd.dist, sample_data(phy97.f.c.coral.Plob)$Dist,  draw = c("polygon"),
#          col = sitecols, alpha=0.2, lwd=0.05)
# legend("topleft", bty="n", pch=1:6,
#        legend=levels(as.factor(sample_data(phy97.f.c.coral.Plob)$Dist)),
#        col=sitecols)
# 
# boxplot(sediment.bd.dist,col=sitecols,ylim=c(0,0.12))
# boxplot(water.bd.dist,col=sitecols,ylim=c(0,0.12))
# boxplot(Peyd.bd.dist,col=sitecols,ylim=c(0,0.12))
# boxplot(MAeq.bd.dist,col=sitecols,ylim=c(0,0.12))
# boxplot(Plob.bd.dist,col=sitecols,ylim=c(0,0.12))
# 
# 
# dev.off()

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
legend("topleft", bty="n", pch=1:6, cex=2.5,
       legend=levels(as.factor(sample_data(phy97.f.c.coral.Peyd)$Dist)),
       col=sitecols)
text(-0.089,-0.037,paste("F = ", round((betadisper.Peyd$`F value`[1]),
                                     digits=1)),cex=2)
text(-0.089,-0.046,paste("p ",(ifelse(betadisper.Peyd$`Pr(>F)`[1]<0.001,"< 0.001", paste("= ",round(betadisper.Peyd$`Pr(>F)`[1],digits=3))))),cex=2)

# Montipora
plot(MAeq.bd.dist, hull=F, label=F, 
     main=expression(paste(italic("Montipora aequituburculata"))), col=sitecols,
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


dev.off()

# 
# 
# jpeg(filename="figures/pcoa_by_dist_boxplots.jpg", 
#      width = 10, height = 6, units="in",res = 300)
# par(mfrow=c(1,5))
# par(mar=c(4.5,4.5,1,1))
# layout(matrix(c(1,2,3,4,5), 1, 5, byrow = TRUE), 
#        widths=c(1.32,1,1,1,1), heights=c(1))
# 
# boxplot(sediment.bd.dist,col=sitecols,ylim=c(0,0.12),axes=FALSE,frame.plot=TRUE,cex.axis=1.5,ylab="Distance to Centroid",cex.lab=1.5)
# axis(side = 1, at=c(1,2), labels = c("Medium","Very High"),cex.axis=1.25)
# axis(side = 1, at=c(1.5), labels = c("Sediment"),line=1.5,tck=FALSE,cex.axis=1.5)
# axis(side = 2, cex.axis=1.5)
# par(mar=c(4.5,0.25,1,1))
# boxplot(water.bd.dist,col=sitecols,ylim=c(0,0.12),axes=FALSE,frame.plot=TRUE,yaxt='',ylab="")
# axis(side = 1, at=c(1,2), labels = c("Medium","Very High"),cex.axis=1.25)
# axis(side = 1, at=c(1.55), labels = c("Water"),line=1.5,tck=FALSE,cex.axis=1.5)
# boxplot(Peyd.bd.dist,col=sitecols,ylim=c(0,0.12),axes=FALSE,frame.plot=TRUE,yaxt='n',ylab="")
# axis(side = 1, at=c(1,2), labels = c("Medium","Very High"),cex.axis=1.25)
# axis(side = 1, at=c(1.5), labels = c("P. eydouxi"),line=1.5,tck=FALSE,cex.axis=1.5)
# boxplot(MAeq.bd.dist,col=sitecols,ylim=c(0,0.12),axes=FALSE,frame.plot=TRUE,yaxt='n',ylab="")
# axis(side = 1, at=c(1,2), labels = c("Medium","Very High"),cex.axis=1.25)
# axis(side = 1, at=c(1.5), labels = c("M. aequituberculata"),line=1.5,tck=FALSE,cex.axis=1.5)
# boxplot(Plob.bd.dist,col=sitecols,ylim=c(0,0.12),axes=FALSE,frame.plot=TRUE,yaxt='n',ylab="")
# axis(side = 1, at=c(1,2), labels = c("Medium","Very High"),cex.axis=1.25)
# axis(side = 1, at=c(1.5), labels = c("P. lobata"),line=1.5,tck=FALSE,cex.axis=1.5)
# 
# 
# dev.off()
# 
# sediment.metadata <- as(sample_data(phy97.f.c.sediment), "data.frame")
# adonis(distance(phy97.f.c.sediment, method="wunifrac") ~ field_season + site,
#        data = sediment.metadata)
# 
# water.metadata <- as(sample_data(phy97.f.c.water), "data.frame")
# adonis(distance(phy97.f.c.water, method="wunifrac") ~ field_season + site,
#        data = water.metadata)
# 
# coral.sp.metadata <- as(sample_data(phy97.f.c.coral), "data.frame")
# adonis(distance(phy97.f.c.coral, method="wunifrac") ~ Coral_Species + field_season + site,
#        data = coral.sp.metadata)
# 
# MAeq.metadata <- as(sample_data(phy97.f.c.coral.MAeq), "data.frame")
# Peyd.metadata <- as(sample_data(phy97.f.c.coral.Peyd), "data.frame")
# Plob.metadata <- as(sample_data(phy97.f.c.coral.Plob), "data.frame")
# 
# 
# adonis(distance(phy97.f.c.sediment, method="wunifrac") ~ field_season + Dist,
#        data = sediment.metadata)
# adonis(distance(phy97.f.c.water, method="wunifrac") ~ field_season + Dist,
#        data = water.metadata)
# adonis(distance(phy97.f.c.coral, method="wunifrac") ~ Coral_Species + field_season + Dist,
#        data = coral.sp.metadata)
# 
# adonis(distance(phy97.f.c.sediment, method="wunifrac") ~ Dist,
#        data = sediment.metadata)
# adonis(distance(phy97.f.c.water, method="wunifrac") ~ Dist,
#        data = water.metadata)
# adonis(distance(phy97.f.c.coral, method="wunifrac") ~ Dist,
#        data = coral.sp.metadata)
# adonis(distance(phy97.f.c.coral.MAeq, method="wunifrac") ~ Dist,
#        data = MAeq.metadata)
# adonis(distance(phy97.f.c.coral.Peyd, method="wunifrac") ~ Dist,
#        data = Peyd.metadata)
# adonis(distance(phy97.f.c.coral.Plob, method="wunifrac") ~ Dist,
#        data = Plob.metadata)
# 
# ord1 <- ordinate(physeq = phy97.f.c.sediment,
#                  method = "CAP",distance = "wunifrac",
#                  formula = ~ Dist + field_season + Dist_fs)
# finalmodel.ord1 <- ordistep(ord1,formula = ~ Dist + field_season + Dist_fs,
#                             direction = c("both"), Pin = 0.05, Pout = 0.1,
#                             pstep=100, perm.max=1000, steps=50, trace=TRUE)
# 
# plot_ordination(physeq = phy97.f.c.sediment,ordination = ord1,shape="Dist",color="Dist_fs")+
#   scale_color_manual(values = c("blue","purple","lightblue","red","maroon","orange"))+ordi
#   ordihull(ord1, sample_data(phy97.f.c.sediment)$field_season,  draw = c("polygon"),
#            col = timecols, alpha=0.2, lwd=0.05)
# 
# ordiplot(ord=ord1,type="points",display="sites")
# ordihull(ord1, sample_data(phy97.f.c.sediment)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b="lightblue",HighMed_KI2014="purple",
#                  HighMed_KI2015a_Post="blue",VeryHigh_KI2015a_Post="darkred",
#                  VeryHigh_KI2015b="maroon",VeryHigh_KI2014="orange"), alpha=0.2, lwd=0.05)
# ordihull(ord1, sample_data(phy97.f.c.sediment)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b="lightblue",HighMed_KI2014="purple",
#                  HighMed_KI2015a_Post="blue",VeryHigh_KI2015a_Post="darkred",
#                  VeryHigh_KI2015b="maroon",VeryHigh_KI2014="orange"), alpha=0.2, lwd=0.05)
# 
# sample_data(phy97.f.c.sediment)$Dist_fs <- paste(sample_data(phy97.f.c.sediment)$Dist,sample_data(phy97.f.c.sediment)$field_season,sep="_")
# 
# ordiplot(ord=ord1,type="points",display="sites")
# ordihull(ord1, sample_data(phy97.f.c.sediment)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b=sitecols[[1]],HighMed_KI2014=sitecols[[1]],
#                  HighMed_KI2015a_Post=sitecols[[1]],VeryHigh_KI2015a_Post=sitecols[[2]],
#                  VeryHigh_KI2015b=sitecols[[2]],VeryHigh_KI2014=sitecols[[2]]), 
#          alpha=0.5, lwd=0.05)
# ordihull(ord1, sample_data(phy97.f.c.sediment)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b=timecols[[4]],HighMed_KI2014=timecols[[1]],
#                  HighMed_KI2015a_Post=timecols[[3]],VeryHigh_KI2015a_Post=timecols[[3]],
#                  VeryHigh_KI2015b=timecols[[4]],VeryHigh_KI2014=timecols[[1]]), alpha=0.2,lwd=3)
# 
# 
# jpeg(filename="figures/pcoa_by_dist_and_time_sediment.jpg", 
#      width = 8, height = 8, units="in",res = 300)
# ordiplot(ord=ord1,type="points",display="sites")
# ordihull(ord1, sample_data(phy97.f.c.sediment)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b=sitecols[[1]],HighMed_KI2014=sitecols[[1]],
#                  HighMed_KI2015a_Post=sitecols[[1]],VeryHigh_KI2015a_Post=sitecols[[2]],
#                  VeryHigh_KI2015b=sitecols[[2]],VeryHigh_KI2014=sitecols[[2]]), 
#          alpha=0.5, lwd=0.05)
# ordihull(ord1, sample_data(phy97.f.c.sediment)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b=timecols[[4]],HighMed_KI2014=timecols[[1]],
#                  HighMed_KI2015a_Post=timecols[[3]],VeryHigh_KI2015a_Post=timecols[[3]],
#                  VeryHigh_KI2015b=timecols[[4]],VeryHigh_KI2014=timecols[[1]]), alpha=0.2,lwd=3)
# dev.off()
# 
# 
# sample_data(phy97.f.c.coral.Peyd)$Dist_fs <- paste(sample_data(phy97.f.c.coral.Peyd)$Dist,sample_data(phy97.f.c.coral.Peyd)$field_season,sep="_")
# 
# sample_sums(phy97.f.c.coral.Peyd)
# 
# ord2 <- ordinate(physeq = phy97.f.c.coral.Peyd,
#                  method = "CAP",distance = "wunifrac",
#                  formula = ~ Dist + field_season + Dist_fs)
# 
# 
# jpeg(filename="figures/pcoa_by_dist_and_time_Peyd.jpg", 
#      width = 8, height = 8, units="in",res = 300)
# ordiplot(ord=ord2,type="points",display="sites")
# ordihull(ord2, sample_data(phy97.f.c.coral.Peyd)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b=sitecols[[1]],HighMed_KI2014=sitecols[[1]],
#                  HighMed_KI2015a_Post=sitecols[[1]],VeryHigh_KI2015a_Post=sitecols[[2]],
#                  VeryHigh_KI2015b=sitecols[[2]],VeryHigh_KI2014=sitecols[[2]]), 
#          alpha=0.5, lwd=0.05)
# ordihull(ord2, sample_data(phy97.f.c.coral.Peyd)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b=timecols[[4]],HighMed_KI2014=timecols[[1]],
#                  HighMed_KI2015a_Post=timecols[[3]],VeryHigh_KI2015a_Post=timecols[[3]],
#                  VeryHigh_KI2015b=timecols[[4]],VeryHigh_KI2014=timecols[[1]]), alpha=0.2,lwd=3)
# dev.off()
# 
# 
# 
# sample_data(phy97.f.c.coral.MAeq)$Dist_fs <- paste(sample_data(phy97.f.c.coral.MAeq)$Dist,sample_data(phy97.f.c.coral.MAeq)$field_season,sep="_")
# 
# ord3 <- ordinate(physeq = phy97.f.c.coral.MAeq,
#                  method = "CAP",distance = "wunifrac",
#                  formula = ~ Dist + field_season + Dist_fs)
# 
# 
# jpeg(filename="figures/pcoa_by_dist_and_time_MAeq.jpg", 
#      width = 8, height = 8, units="in",res = 300)
# ordiplot(ord=ord3,type="points",display="sites")
# ordihull(ord3, sample_data(phy97.f.c.coral.MAeq)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b=sitecols[[1]],HighMed_KI2014=sitecols[[1]],
#                  HighMed_KI2015a_Post=sitecols[[1]],VeryHigh_KI2015a_Post=sitecols[[2]],
#                  VeryHigh_KI2015b=sitecols[[2]],VeryHigh_KI2014=sitecols[[2]]), 
#          alpha=0.5, lwd=0.05)
# ordihull(ord3, sample_data(phy97.f.c.coral.MAeq)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b=timecols[[4]],HighMed_KI2014=timecols[[1]],
#                  HighMed_KI2015a_Post=timecols[[3]],VeryHigh_KI2015a_Post=timecols[[3]],
#                  VeryHigh_KI2015b=timecols[[4]],VeryHigh_KI2014=timecols[[1]]), alpha=0.2,lwd=3)
# dev.off()
# 
# 
# sample_data(phy97.f.c.water)$Dist_fs <- paste(sample_data(phy97.f.c.water)$Dist,sample_data(phy97.f.c.water)$field_season,sep="_")
# 
# ord5 <- ordinate(physeq = phy97.f.c.water,
#                  method = "CAP",distance = "wunifrac",
#                  formula = ~ Dist + field_season + Dist_fs)
# 
# 
# jpeg(filename="figures/pcoa_by_dist_and_time_water.jpg", 
#      width = 8, height = 8, units="in",res = 300)
# ordiplot(ord=ord5,type="points",display="sites")
# ordihull(ord5, sample_data(phy97.f.c.water)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b=sitecols[[1]],HighMed_KI2014=sitecols[[1]],
#                  HighMed_KI2015a_Post=sitecols[[1]],VeryHigh_KI2015a_Post=sitecols[[2]],
#                  VeryHigh_KI2015b=sitecols[[2]],VeryHigh_KI2014=sitecols[[2]]), 
#          alpha=0.5, lwd=0.05)
# ordihull(ord5, sample_data(phy97.f.c.water)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b=timecols[[4]],HighMed_KI2014=timecols[[1]],
#                  HighMed_KI2015a_Post=timecols[[3]],VeryHigh_KI2015a_Post=timecols[[3]],
#                  VeryHigh_KI2015b=timecols[[4]],VeryHigh_KI2014=timecols[[1]]), alpha=0.2,lwd=3)
# dev.off()
# 
# 
# jpeg(filename="figures/pcoa_by_dist_and_time_Plob.jpg", 
#      width = 8, height = 8, units="in",res = 300)
# ordiplot(ord=ord4,type="points",display="sites")
# ordihull(ord4, sample_data(phy97.f.c.water)$Dist_fs,  draw = c("polygon"),
#          col = c(HighMed_KI2015b=sitecols[[1]],HighMed_KI2014=sitecols[[1]],
#                  HighMed_KI2015a_Post=sitecols[[1]],VeryHigh_KI2015a_Post=sitecols[[2]],
#                  VeryHigh_KI2015b=sitecols[[2]],VeryHigh_KI2014=sitecols[[2]]), 
#          alpha=0.5, lwd=0.05)
# ordihull(ord4, sample_data(phy97.f.c.water)$Dist_fs,  draw = c("lines"),
#          col = c(HighMed_KI2015b=timecols[[4]],HighMed_KI2014=timecols[[1]],
#                  HighMed_KI2015a_Post=timecols[[3]],VeryHigh_KI2015a_Post=timecols[[3]],
#                  VeryHigh_KI2015b=timecols[[4]],VeryHigh_KI2014=timecols[[1]]), alpha=0.2,lwd=3)
# dev.off()
