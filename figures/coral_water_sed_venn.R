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


phy97.f.c.coral.VH <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Dist=="VeryHigh")
phy97.f.c.coral.M <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Dist=="HighMed")
phy97.f.c.sediment.VH <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$Dist=="VeryHigh")
phy97.f.c.sediment.M <- subset_samples(phy97.f.c.sediment,sample_data(phy97.f.c.sediment)$Dist=="HighMed")
phy97.f.c.water.VH <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$Dist=="VeryHigh")
phy97.f.c.water.M <- subset_samples(phy97.f.c.water,sample_data(phy97.f.c.water)$Dist=="HighMed")

sediment.VH.types <- unique(data.frame(tax_table(phy97.f.c.sediment.VH))$hit)
sediment.M.types <- unique(data.frame(tax_table(phy97.f.c.sediment.M))$hit)
water.VH.types <- unique(data.frame(tax_table(phy97.f.c.water.VH))$hit)
water.M.types <- unique(data.frame(tax_table(phy97.f.c.water.M))$hit)
coral.VH.types <- unique(data.frame(tax_table(phy97.f.c.coral.VH))$hit)
coral.M.types <- unique(data.frame(tax_table(phy97.f.c.coral.M))$hit)

sediment.VH.types.subclade <- sediment.VH.types
sediment.VH.types.subclade <- gsub("_.*","",sediment.VH.types.subclade)
sediment.VH.types.subclade <- gsub("\\..*","",sediment.VH.types.subclade)
sediment.VH.types.subclade <- unique(sediment.VH.types.subclade)
sediment.M.types.subclade <- sediment.M.types
sediment.M.types.subclade <- gsub("_.*","",sediment.M.types.subclade)
sediment.M.types.subclade <- gsub("\\..*","",sediment.M.types.subclade)
sediment.M.types.subclade <- unique(sediment.M.types.subclade)

water.VH.types.subclade <- water.VH.types
water.VH.types.subclade <- gsub("_.*","",water.VH.types.subclade)
water.VH.types.subclade <- gsub("\\..*","",water.VH.types.subclade)
water.VH.types.subclade <- unique(water.VH.types.subclade)
water.M.types.subclade <- water.M.types
water.M.types.subclade <- gsub("_.*","",water.M.types.subclade)
water.M.types.subclade <- gsub("\\..*","",water.M.types.subclade)
water.M.types.subclade <- unique(water.M.types.subclade)

coral.VH.types.subclade <- coral.VH.types
coral.VH.types.subclade <- gsub("_.*","",coral.VH.types.subclade)
coral.VH.types.subclade <- gsub("\\..*","",coral.VH.types.subclade)
coral.VH.types.subclade <- unique(coral.VH.types.subclade)
coral.M.types.subclade <- coral.M.types
coral.M.types.subclade <- gsub("_.*","",coral.M.types.subclade)
coral.M.types.subclade <- gsub("\\..*","",coral.M.types.subclade)
coral.M.types.subclade <- unique(coral.M.types.subclade)

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
VennDiag3 <- euler(c("Coral" = (c.VH), "Sediment" = (s.VH), "Water" = (w.VH), "Coral&Sediment" = (cs.only.VH), "Coral&Water" = (cw.only.VH), "Sediment&Water" = (sw.only.VH), "Coral&Sediment&Water" = (csw.length.VH)))
plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
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
VennDiag3 <- euler(c("Coral" = (c.M), "Sediment" = (s.M), "Water" = (w.M), "Coral&Sediment" = (cs.only.M), "Coral&Water" = (cw.only.M), "Sediment&Water" = (sw.only.M), "Coral&Sediment&Water" = (csw.length.M)))
plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()
