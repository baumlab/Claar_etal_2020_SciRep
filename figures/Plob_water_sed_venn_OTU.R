rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

compcols["coral"] <- speccols["P. lobata"]

cs <- intersect(Plob.denovo.subclade,sediment.denovo.subclade)
csw <- intersect(cs,water.denovo.subclade)
cs.only <- length(cs)-length(csw)
csw.length <- length(csw)
cw <- intersect(Plob.denovo.subclade,water.denovo.subclade)
cw.only <- length(cw)-length(csw)
sw <- intersect(sediment.denovo.subclade,water.denovo.subclade)
sw.only <- length(sw)-length(csw)
c <- length(Plob.denovo.subclade) - cs.only - cw.only - csw.length
s <- length(sediment.denovo.subclade) - sw.only - cs.only - csw.length
w <- length(water.denovo.subclade) - sw.only - cw.only - csw.length

jpeg(filename="figures/Plob_water_sed_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag7 <- euler(c("Plob" = (c), "Sediment" = (s), "Water" = (w), "Plob&Sediment" = (cs.only), "Plob&Water" = (cw.only), "Sediment&Water" = (sw.only), "Plob&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

pdf(file="figures/Plob_water_sed_venn_OTU.pdf", 
    width = 4, height = 4)
VennDiag7 <- euler(c("Plob" = (c), "Sediment" = (s), "Water" = (w), "Plob&Sediment" = (cs.only), "Plob&Water" = (cw.only), "Sediment&Water" = (sw.only), "Plob&Sediment&Water" = (csw.length)))
plot(VennDiag7, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()
