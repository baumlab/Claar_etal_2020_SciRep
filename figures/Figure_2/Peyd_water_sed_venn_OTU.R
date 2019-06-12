rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

compcols["coral"] <- speccols["P. eydouxi"]

cs <- intersect(Peyd.ASVs,sediment.ASVs)
csw <- intersect(cs,water.ASVs)
cs.only <- length(cs)-length(csw)
csw.length <- length(csw)
cw <- intersect(Peyd.ASVs,water.ASVs)
cw.only <- length(cw)-length(csw)
sw <- intersect(sediment.ASVs,water.ASVs)
sw.only <- length(sw)-length(csw)
c <- length(Peyd.ASVs) - cs.only - cw.only - csw.length
s <- length(sediment.ASVs) - sw.only - cs.only - csw.length
w <- length(water.ASVs) - sw.only - cw.only - csw.length

jpeg(filename="figures/Figure_2/Peyd_water_sed_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag7 <- euler(c("Peyd" = (c), "Sediment" = (s), "Water" = (w), "Peyd&Sediment" = (cs.only), "Peyd&Water" = (cw.only), "Sediment&Water" = (sw.only), "Peyd&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

pdf(file="figures/Figure_2/Peyd_water_sed_venn_OTU.pdf", 
    width = 4, height = 4)
VennDiag7 <- euler(c("Peyd" = (c), "Sediment" = (s), "Water" = (w), "Peyd&Sediment" = (cs.only), "Peyd&Water" = (cw.only), "Sediment&Water" = (sw.only), "Peyd&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()
