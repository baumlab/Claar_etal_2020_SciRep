rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

cs <- intersect(coral.types.subclade,sediment.types.subclade)
csw <- intersect(cs,water.types.subclade)
cs.only <- length(cs)-length(csw)
csw.length <- length(csw)
cw <- intersect(coral.types.subclade,water.types.subclade)
cw.only <- length(cw)-length(csw)
sw <- intersect(sediment.types.subclade,water.types.subclade)
sw.only <- length(sw)-length(csw)
c <- length(coral.types.subclade) - cs.only - cw.only - csw.length
s <- length(sediment.types.subclade) - sw.only - cs.only - csw.length
w <- length(water.types.subclade) - sw.only - cw.only - csw.length

jpeg(filename="figures/coral_water_sed_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag3 <- euler(c("Coral" = (c), "Sediment" = (s), "Water" = (w), "Coral&Sediment" = (cs.only), "Coral&Water" = (cw.only), "Sediment&Water" = (sw.only), "Coral&Sediment&Water" = (csw.length)))
plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

write.table(file="analyses/coral_subclades.csv",coral.types.subclade,row.names = F,col.names = F,sep=",",quote=F)
write.table(file="analyses/sediment_subclades.csv",sediment.types.subclade,row.names = F,col.names = F,sep=",",quote=F)
write.table(file="analyses/water_subclades.csv",water.types.subclade,row.names = F,col.names = F,sep=",",quote=F)

coral.sediment.subclades <- coral.types.subclade[which((coral.types.subclade %in% sediment.types.subclade) & !(coral.types.subclade %in% water.types.subclade))]
water.sediment.subclades <- water.types.subclade[which((water.types.subclade %in% sediment.types.subclade) & !(water.types.subclade %in% coral.types.subclade))]
coral.water.subclades <- coral.types.subclade[which((coral.types.subclade %in% water.types.subclade) & !(coral.types.subclade %in% sediment.types.subclade))]
coral.only.subclades <- coral.types.subclade[which(!(coral.types.subclade %in% sediment.types.subclade) & !(coral.types.subclade %in% water.types.subclade))]
water.only.subclades <- water.types.subclade[which(!(water.types.subclade %in% sediment.types.subclade) & !(water.types.subclade %in% coral.types.subclade))]
sediment.only.subclades <- sediment.types.subclade[which(!(sediment.types.subclade %in% coral.types.subclade) & !(sediment.types.subclade %in% water.types.subclade))]
coral.sediment.water.subclades <- coral.types.subclade[which((coral.types.subclade %in% sediment.types.subclade) & (coral.types.subclade %in% water.types.subclade))]


cs.VH <- intersect(coral.VH.types.subclade,sediment.VH.types.subclade)
csw.VH <- intersect(cs.VH,water.VH.types.subclade)
cs.only.VH <- length(cs.VH)-length(csw.VH)
csw.length.VH <- length(csw.VH)
cw.VH <- intersect(coral.VH.types.subclade,water.VH.types.subclade)
cw.only.VH <- length(cw.VH)-length(csw.VH)
sw.VH <- intersect(sediment.VH.types.subclade,water.VH.types.subclade)
sw.only.VH <- length(sw.VH)-length(csw.VH)
c.VH <- length(coral.VH.types.subclade) - cs.only.VH - cw.only.VH - csw.length.VH
s.VH <- length(sediment.VH.types.subclade) - sw.only.VH - cs.only.VH - csw.length.VH
w.VH <- length(water.VH.types.subclade) - sw.only.VH - cw.only.VH - csw.length.VH

jpeg(filename="figures/coral_water_sed_venn_VH.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag5 <- euler(c("Coral" = (c.VH), "Sediment" = (s.VH), "Water" = (w.VH), "Coral&Sediment" = (cs.only.VH), "Coral&Water" = (cw.only.VH), "Sediment&Water" = (sw.only.VH), "Coral&Sediment&Water" = (csw.length.VH)))
plot(VennDiag5, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

cs.M <- intersect(coral.M.types.subclade,sediment.M.types.subclade)
csw.M <- intersect(cs.M,water.M.types.subclade)
cs.only.M <- length(cs.M)-length(csw.M)
csw.length.M <- length(csw.M)
cw.M <- intersect(coral.M.types.subclade,water.M.types.subclade)
cw.only.M <- length(cw.M)-length(csw.M)
sw.M <- intersect(sediment.M.types.subclade,water.M.types.subclade)
sw.only.M <- length(sw.M)-length(csw.M)
c.M <- length(coral.M.types.subclade) - cs.only.M - cw.only.M - csw.length.M
s.M <- length(sediment.M.types.subclade) - sw.only.M - cs.only.M - csw.length.M
w.M <- length(water.M.types.subclade) - sw.only.M - cw.only.M - csw.length.M

jpeg(filename="figures/coral_water_sed_venn_M.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag6 <- euler(c("Coral" = (c.M), "Sediment" = (s.M), "Water" = (w.M), "Coral&Sediment" = (cs.only.M), "Coral&Water" = (cw.only.M), "Sediment&Water" = (sw.only.M), "Coral&Sediment&Water" = (csw.length.M)))
plot(VennDiag6, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

vd3 <- plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
vd5 <- plot(VennDiag5, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
vd6 <- plot(VennDiag6, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))


jpeg(filename="figures/coral_water_sed_venn_ALL_VH_M.jpg", 
     width = 12, height = 4, units="in",res = 300)
grid.arrange(vd3,vd6,vd5,ncol=3)
dev.off()
