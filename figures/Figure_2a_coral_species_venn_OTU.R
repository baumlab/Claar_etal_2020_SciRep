rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

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


jpeg(filename="figures/Figure_2a_coral_species_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag4 <- euler(c("P. grandis" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. grandis&P. lobata" = (PeydPlob.only), "P. grandis&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off()