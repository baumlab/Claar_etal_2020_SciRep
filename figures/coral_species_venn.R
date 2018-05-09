rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")


PeydPlob <- intersect(Peyd.denovo.subclade,Peyd.denovo.subclade)
PeydPlobMaeq <- intersect(PeydPlob,Peyd.denovo.subclade)
PeydPlob.only <- length(PeydPlob)-length(PeydPlobMaeq)
PeydPlobMaeq.length <- length(PeydPlobMaeq)
PeydMaeq <- intersect(Peyd.denovo.subclade,Peyd.denovo.subclade)
PeydMaeq.only <- length(PeydMaeq)-length(PeydPlobMaeq)
PlobMaeq <- intersect(Peyd.denovo.subclade,Peyd.denovo.subclade)
PlobMaeq.only <- length(PlobMaeq)-length(PeydPlobMaeq)
Peyd <- length(Peyd.denovo.subclade) - PeydPlob.only - PeydMaeq.only - PeydPlobMaeq.length
Plob <- length(Peyd.denovo.subclade) - PlobMaeq.only - PeydPlob.only - PeydPlobMaeq.length
Maeq <- length(Peyd.denovo.subclade) - PlobMaeq.only - PeydMaeq.only - PeydPlobMaeq.length


jpeg(filename="figures/coral_species_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag4 <- euler(c("P. grandis" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. grandis&P. lobata" = (PeydPlob.only), "P. grandis&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
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
VennDiag4 <- euler(c("P. grandis" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. grandis&P. lobata" = (PeydPlob.only), "P. grandis&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off()


PeydPlob.VH <- intersect(Peyd.VH.denovo.subclade,Plob.VH.denovo.subclade)
PeydPlobMaeq.VH <- intersect(PeydPlob.VH,MAeq.VH.denovo.subclade)
PeydPlob.only.VH <- length(PeydPlob.VH)-length(PeydPlobMaeq.VH)
PeydPlobMaeq.length.VH <- length(PeydPlobMaeq.VH)
PeydMaeq.VH <- intersect(Peyd.VH.denovo.subclade,MAeq.VH.denovo.subclade)
PeydMaeq.only.VH <- length(PeydMaeq.VH)-length(PeydPlobMaeq.VH)
PlobMaeq.VH <- intersect(Plob.VH.denovo.subclade,MAeq.VH.denovo.subclade)
PlobMaeq.only.VH <- length(PlobMaeq.VH)-length(PeydPlobMaeq.VH)
Peyd.VH <- length(Peyd.VH.denovo.subclade) - PeydPlob.only.VH - PeydMaeq.only.VH - PeydPlobMaeq.length.VH
Plob.VH <- length(Plob.VH.denovo.subclade) - PlobMaeq.only.VH - PeydPlob.only.VH - PeydPlobMaeq.length.VH
Maeq.VH <- length(MAeq.VH.denovo.subclade) - PlobMaeq.only.VH - PeydMaeq.only.VH - PeydPlobMaeq.length.VH


jpeg(filename="figures/coral_species_venn_OTU_VH.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag10 <- euler(c("P. grandis" = (Peyd.VH), "P. lobata" = (Plob.VH), "M. aequituberculata" = (Maeq.VH), "P. grandis&P. lobata" = (PeydPlob.only.VH), "P. grandis&M. aequituberculata" = (PeydMaeq.only.VH), "P. lobata&M. aequituberculata" = (PlobMaeq.only.VH), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length.VH)))
plot(VennDiag10, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off()


PeydPlob.M <- intersect(Peyd.M.denovo.subclade,Plob.M.denovo.subclade)
PeydPlobMaeq.M <- intersect(PeydPlob.M,MAeq.M.denovo.subclade)
PeydPlob.only.M <- length(PeydPlob.M)-length(PeydPlobMaeq.M)
PeydPlobMaeq.length.M <- length(PeydPlobMaeq.M)
PeydMaeq.M <- intersect(Peyd.M.denovo.subclade,MAeq.M.denovo.subclade)
PeydMaeq.only.M <- length(PeydMaeq.M)-length(PeydPlobMaeq.M)
PlobMaeq.M <- intersect(Plob.M.denovo.subclade,MAeq.M.denovo.subclade)
PlobMaeq.only.M <- length(PlobMaeq.M)-length(PeydPlobMaeq.M)
Peyd.M <- length(Peyd.M.denovo.subclade) - PeydPlob.only.M - PeydMaeq.only.M - PeydPlobMaeq.length.M
Plob.M <- length(Plob.M.denovo.subclade) - PlobMaeq.only.M - PeydPlob.only.M - PeydPlobMaeq.length.M
Maeq.M <- length(MAeq.M.denovo.subclade) - PlobMaeq.only.M - PeydMaeq.only.M - PeydPlobMaeq.length.M


jpeg(filename="figures/coral_species_venn_OTU_M.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag11 <- euler(c("P. grandis" = (Peyd.M), "P. lobata" = (Plob.M), "M. aequituberculata" = (Maeq.M), "P. grandis&P. lobata" = (PeydPlob.only.M), "P. grandis&M. aequituberculata" = (PeydMaeq.only.M), "P. lobata&M. aequituberculata" = (PlobMaeq.only.M), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length.M)))
plot(VennDiag11, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off()


vd4 <- plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
vd10 <- plot(VennDiag10, quantities = TRUE, font=1, cex=1, alpha=0.5,
             fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
vd11 <- plot(VennDiag11, quantities = TRUE, font=1, cex=1, alpha=0.5,
             fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))

jpeg(filename="figures/coral_species_venn_OTU_ALL_VH_M.jpg", 
     width = 12, height = 4, units="in",res = 300)
grid.arrange(vd4,vd10,vd11,ncol=3)
dev.off()
