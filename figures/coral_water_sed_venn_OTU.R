rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")


cs <- intersect(coral.denovo.subclade,sediment.denovo.subclade)
csw <- intersect(cs,water.denovo.subclade)
cs.only <- length(cs)-length(csw)
csw.length <- length(csw)
cw <- intersect(coral.denovo.subclade,water.denovo.subclade)
cw.only <- length(cw)-length(csw)
sw <- intersect(sediment.denovo.subclade,water.denovo.subclade)
sw.only <- length(sw)-length(csw)
c <- length(coral.denovo.subclade) - cs.only - cw.only - csw.length
s <- length(sediment.denovo.subclade) - sw.only - cs.only - csw.length
w <- length(water.denovo.subclade) - sw.only - cw.only - csw.length

jpeg(filename="figures/coral_water_sed_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag7 <- euler(c("Coral" = (c), "Sediment" = (s), "Water" = (w), "Coral&Sediment" = (cs.only), "Coral&Water" = (cw.only), "Sediment&Water" = (sw.only), "Coral&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()


PeydPlob <- intersect(Peyd.denovo.subclade,Plob.denovo.subclade)
PeydPlobMaeq <- intersect(PeydPlob,MAeq.denovo.subclade)
PeydPlob.only <- length(PeydPlob)-length(PeydPlobMaeq)
PeydPlobMaeq.length <- length(PeydPlobMaeq)
PeydMaeq <- intersect(Peyd.denovo.subclade,MAeq.denovo.subclade)
PeydMaeq.only <- length(PeydMaeq)-length(PeydPlobMaeq)
PlobMaeq <- intersect(Plob.denovo.subclade,MAeq.denovo.subclade)
PlobMaeq.only <- length(PlobMaeq)-length(PeydPlobMaeq)
Peyd <- length(Peyd.denovo.subclade) - PeydPlob.only - PeydMaeq.only - PeydPlobMaeq.length
Plob <- length(Plob.denovo.subclade) - PlobMaeq.only - PeydPlob.only - PeydPlobMaeq.length
Maeq <- length(MAeq.denovo.subclade) - PlobMaeq.only - PeydMaeq.only - PeydPlobMaeq.length


jpeg(filename="figures/coral_species_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag4 <- euler(c("P. eydouxi" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. eydouxi&P. lobata" = (PeydPlob.only), "P. eydouxi&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. eydouxi&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off()


cs.VH <- intersect(coral.VH.denovo.subclade,sediment.VH.denovo.subclade)
csw.VH <- intersect(cs.VH,water.VH.denovo.subclade)
cs.only.VH <- length(cs.VH)-length(csw.VH)
csw.length.VH <- length(csw.VH)
cw.VH <- intersect(coral.VH.denovo.subclade,water.VH.denovo.subclade)
cw.only.VH <- length(cw.VH)-length(csw.VH)
sw.VH <- intersect(sediment.VH.denovo.subclade,water.VH.denovo.subclade)
sw.only.VH <- length(sw.VH)-length(csw.VH)
c.VH <- length(coral.VH.denovo.subclade) - cs.only.VH - cw.only.VH - csw.length.VH
s.VH <- length(sediment.VH.denovo.subclade) - sw.only.VH - cs.only.VH - csw.length.VH
w.VH <- length(water.VH.denovo.subclade) - sw.only.VH - cw.only.VH - csw.length.VH

jpeg(filename="figures/coral_water_sed_venn_VH_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag5 <- euler(c("Coral" = (c.VH), "Sediment" = (s.VH), "Water" = (w.VH), "Coral&Sediment" = (cs.only.VH), "Coral&Water" = (cw.only.VH), "Sediment&Water" = (sw.only.VH), "Coral&Sediment&Water" = (csw.length.VH)))
plot(VennDiag5, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

cs.M <- intersect(coral.M.denovo.subclade,sediment.M.denovo.subclade)
csw.M <- intersect(cs.M,water.M.denovo.subclade)
cs.only.M <- length(cs.M)-length(csw.M)
csw.length.M <- length(csw.M)
cw.M <- intersect(coral.M.denovo.subclade,water.M.denovo.subclade)
cw.only.M <- length(cw.M)-length(csw.M)
sw.M <- intersect(sediment.M.denovo.subclade,water.M.denovo.subclade)
sw.only.M <- length(sw.M)-length(csw.M)
c.M <- length(coral.M.denovo.subclade) - cs.only.M - cw.only.M - csw.length.M
s.M <- length(sediment.M.denovo.subclade) - sw.only.M - cs.only.M - csw.length.M
w.M <- length(water.M.denovo.subclade) - sw.only.M - cw.only.M - csw.length.M

jpeg(filename="figures/coral_water_sed_venn_M_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag6 <- euler(c("Coral" = (c.M), "Sediment" = (s.M), "Water" = (w.M), "Coral&Sediment" = (cs.only.M), "Coral&Water" = (cw.only.M), "Sediment&Water" = (sw.only.M), "Coral&Sediment&Water" = (csw.length.M)))
plot(VennDiag6, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

vd7 <- plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
vd5 <- plot(VennDiag5, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
vd6 <- plot(VennDiag6, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))


jpeg(filename="figures/coral_water_sed_venn_ALL_VH_M_OTU.jpg", 
     width = 12, height = 4, units="in",res = 300)
grid.arrange(vd7,vd6,vd5,ncol=3)
dev.off()
