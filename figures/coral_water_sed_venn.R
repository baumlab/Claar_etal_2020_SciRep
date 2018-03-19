rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

all.types <- unique(data.frame(tax_table(phy97.f.c))$hit)
sediment.types <- unique(data.frame(tax_table(phy97.f.c.sediment))$hit)
water.types <- unique(data.frame(tax_table(phy97.f.c.water))$hit)
coral.types <- unique(data.frame(tax_table(phy97.f.c.coral))$hit)

sediment.types.subclade <- sediment.types
sediment.types.subclade <- gsub("_.*","",sediment.types.subclade)
sediment.types.subclade <- gsub("\\..*","",sediment.types.subclade)
sediment.types.subclade <- unique(sediment.types.subclade)

water.types.subclade <- water.types
water.types.subclade <- gsub("_.*","",water.types.subclade)
water.types.subclade <- gsub("\\..*","",water.types.subclade)
water.types.subclade <- unique(water.types.subclade)

coral.types.subclade <- coral.types
coral.types.subclade <- gsub("_.*","",coral.types.subclade)
coral.types.subclade <- gsub("\\..*","",coral.types.subclade)
coral.types.subclade <- unique(coral.types.subclade)

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




phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

Peyd.types <- unique(data.frame(tax_table(phy97.f.c.coral.Peyd))$hit)
Peyd.types.subclade <- Peyd.types
Peyd.types.subclade <- gsub("_.*","",Peyd.types.subclade)
Peyd.types.subclade <- gsub("\\..*","",Peyd.types.subclade)
Peyd.types.subclade <- unique(Peyd.types.subclade)

Plob.types <- unique(data.frame(tax_table(phy97.f.c.coral.Plob))$hit)
Plob.types.subclade <- Plob.types
Plob.types.subclade <- gsub("_.*","",Plob.types.subclade)
Plob.types.subclade <- gsub("\\..*","",Plob.types.subclade)
Plob.types.subclade <- unique(Plob.types.subclade)

MAeq.types <- unique(data.frame(tax_table(phy97.f.c.coral.MAeq))$hit)
MAeq.types.subclade <- MAeq.types
MAeq.types.subclade <- gsub("_.*","",MAeq.types.subclade)
MAeq.types.subclade <- gsub("\\..*","",MAeq.types.subclade)
MAeq.types.subclade <- unique(MAeq.types.subclade)

PeydPlob <- intersect(Peyd.types.subclade,Plob.types.subclade)
PeydPlobMaeq <- intersect(PeydPlob,MAeq.types.subclade)
PeydPlob.only <- length(PeydPlob)-length(PeydPlobMaeq)
PeydPlobMaeq.length <- length(PeydPlobMaeq)
PeydMaeq <- intersect(Peyd.types.subclade,MAeq.types.subclade)
PeydMaeq.only <- length(PeydMaeq)-length(PeydPlobMaeq)
PlobMaeq <- intersect(Plob.types.subclade,MAeq.types.subclade)
PlobMaeq.only <- length(PlobMaeq)-length(PeydPlobMaeq)
Peyd <- length(Peyd.types.subclade) - PeydPlob.only - PeydMaeq.only - PeydPlobMaeq.length
Plob <- length(Plob.types.subclade) - PlobMaeq.only - PeydPlob.only - PeydPlobMaeq.length
Maeq <- length(MAeq.types.subclade) - PlobMaeq.only - PeydMaeq.only - PeydPlobMaeq.length


jpeg(filename="figures/coral_species_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag3 <- euler(c("P. eydouxi" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. eydouxi&P. lobata" = (PeydPlob.only), "P. eydouxi&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. eydouxi&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()


