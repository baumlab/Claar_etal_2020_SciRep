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

jpeg(filename="figures/Figure_2b_coral_water_sed_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag7 <- euler(c("Coral" = (c), "Sediment" = (s), "Water" = (w), "Coral&Sediment" = (cs.only), "Coral&Water" = (cw.only), "Sediment&Water" = (sw.only), "Coral&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

pdf(file="figures/Figure_2b_coral_water_sed_venn_OTU.pdf", 
     width = 4, height = 4)
VennDiag7 <- euler(c("Coral" = (c), "Sediment" = (s), "Water" = (w), "Coral&Sediment" = (cs.only), "Coral&Water" = (cw.only), "Sediment&Water" = (sw.only), "Coral&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()
