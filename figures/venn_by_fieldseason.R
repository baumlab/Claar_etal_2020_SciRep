rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")


phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")


sediment.before <- subset_samples(phy97.f.c.sediment, data.frame(sample_data(phy97.f.c.sediment))$field_season == "KI2014",prune=TRUE)
sediment.before <- subset_taxa(sediment.before, taxa_sums(sediment.before) > 0, prune=TRUE)

sediment.storm <- subset_samples(phy97.f.c.sediment, data.frame(sample_data(phy97.f.c.sediment))$field_season == "KI2015a_Post", prune=TRUE)
sediment.storm <- subset_taxa(sediment.storm, taxa_sums(sediment.storm) > 0, prune=TRUE)

sediment.after <- subset_samples(phy97.f.c.sediment, data.frame(sample_data(phy97.f.c.sediment))$field_season == "KI2015b", prune=TRUE)
sediment.after <- subset_taxa(sediment.after, taxa_sums(sediment.after) > 0, prune=TRUE)

sediment.before.types <- unique(data.frame(tax_table(sediment.before))$hit)
sediment.storm.types <- unique(data.frame(tax_table(sediment.storm))$hit)
sediment.after.types <- unique(data.frame(tax_table(sediment.after))$hit)

sediment.before.types.subclade <- sediment.before.types
sediment.before.types.subclade <- gsub("_.*","",sediment.before.types.subclade)
sediment.before.types.subclade <- gsub("\\..*","",sediment.before.types.subclade)
sediment.before.types.subclade <- unique(sediment.before.types.subclade)

sediment.storm.types.subclade <- sediment.storm.types
sediment.storm.types.subclade <- gsub("_.*","",sediment.storm.types.subclade)
sediment.storm.types.subclade <- gsub("\\..*","",sediment.storm.types.subclade)
sediment.storm.types.subclade <- unique(sediment.storm.types.subclade)

sediment.after.types.subclade <- sediment.after.types
sediment.after.types.subclade <- gsub("_.*","",sediment.after.types.subclade)
sediment.after.types.subclade <- gsub("\\..*","",sediment.after.types.subclade)
sediment.after.types.subclade <- unique(sediment.after.types.subclade)

water.before <- subset_samples(phy97.f.c.water, data.frame(sample_data(phy97.f.c.water))$field_season == "KI2014",prune=TRUE)
water.before <- subset_taxa(water.before, taxa_sums(water.before) > 0, prune=TRUE)

water.storm <- subset_samples(phy97.f.c.water, data.frame(sample_data(phy97.f.c.water))$field_season == "KI2015a_Post", prune=TRUE)
water.storm <- subset_taxa(water.storm, taxa_sums(water.storm) > 0, prune=TRUE)

water.after <- subset_samples(phy97.f.c.water, data.frame(sample_data(phy97.f.c.water))$field_season == "KI2015b", prune=TRUE)
water.after <- subset_taxa(water.after, taxa_sums(water.after) > 0, prune=TRUE)

water.before.types <- unique(data.frame(tax_table(water.before))$hit)
water.storm.types <- unique(data.frame(tax_table(water.storm))$hit)
water.after.types <- unique(data.frame(tax_table(water.after))$hit)

water.before.types.subclade <- water.before.types
water.before.types.subclade <- gsub("_.*","",water.before.types.subclade)
water.before.types.subclade <- gsub("\\..*","",water.before.types.subclade)
water.before.types.subclade <- unique(water.before.types.subclade)

water.storm.types.subclade <- water.storm.types
water.storm.types.subclade <- gsub("_.*","",water.storm.types.subclade)
water.storm.types.subclade <- gsub("\\..*","",water.storm.types.subclade)
water.storm.types.subclade <- unique(water.storm.types.subclade)

water.after.types.subclade <- water.after.types
water.after.types.subclade <- gsub("_.*","",water.after.types.subclade)
water.after.types.subclade <- gsub("\\..*","",water.after.types.subclade)
water.after.types.subclade <- unique(water.after.types.subclade)

coral.before <- subset_samples(phy97.f.c.coral, data.frame(sample_data(phy97.f.c.coral))$field_season == "KI2014",prune=TRUE)
coral.before <- subset_taxa(coral.before, taxa_sums(coral.before) > 0, prune=TRUE)

coral.storm <- subset_samples(phy97.f.c.coral, data.frame(sample_data(phy97.f.c.coral))$field_season == "KI2015a_Post", prune=TRUE)
coral.storm <- subset_taxa(coral.storm, taxa_sums(coral.storm) > 0, prune=TRUE)

coral.after <- subset_samples(phy97.f.c.coral, data.frame(sample_data(phy97.f.c.coral))$field_season == "KI2015b", prune=TRUE)
coral.after <- subset_taxa(coral.after, taxa_sums(coral.after) > 0, prune=TRUE)

coral.before.types <- unique(data.frame(tax_table(coral.before))$hit)
coral.storm.types <- unique(data.frame(tax_table(coral.storm))$hit)
coral.after.types <- unique(data.frame(tax_table(coral.after))$hit)

coral.before.types.subclade <- coral.before.types
coral.before.types.subclade <- gsub("_.*","",coral.before.types.subclade)
coral.before.types.subclade <- gsub("\\..*","",coral.before.types.subclade)
coral.before.types.subclade <- unique(coral.before.types.subclade)

coral.storm.types.subclade <- coral.storm.types
coral.storm.types.subclade <- gsub("_.*","",coral.storm.types.subclade)
coral.storm.types.subclade <- gsub("\\..*","",coral.storm.types.subclade)
coral.storm.types.subclade <- unique(coral.storm.types.subclade)

coral.after.types.subclade <- coral.after.types
coral.after.types.subclade <- gsub("_.*","",coral.after.types.subclade)
coral.after.types.subclade <- gsub("\\..*","",coral.after.types.subclade)
coral.after.types.subclade <- unique(coral.after.types.subclade)

# Peyd.before <- subset_samples(phy97.f.c.coral.Peyd, data.frame(sample_data(phy97.f.c.coral.Peyd))$field_season == "KI2014",prune=TRUE)
# Peyd.before <- subset_taxa(Peyd.before, taxa_sums(Peyd.before) > 0, prune=TRUE)
# 
# Peyd.storm <- subset_samples(phy97.f.c.coral.Peyd, data.frame(sample_data(phy97.f.c.coral.Peyd))$field_season == "KI2015a_Post", prune=TRUE)
# Peyd.storm <- subset_taxa(Peyd.storm, taxa_sums(Peyd.storm) > 0, prune=TRUE)
# 
# Peyd.after <- subset_samples(phy97.f.c.coral.Peyd, data.frame(sample_data(phy97.f.c.coral.Peyd))$field_season == "KI2015b", prune=TRUE)
# Peyd.after <- subset_taxa(Peyd.after, taxa_sums(Peyd.after) > 0, prune=TRUE)
# 
# Peyd.before.types <- unique(data.frame(tax_table(Peyd.before))$hit)
# Peyd.storm.types <- unique(data.frame(tax_table(Peyd.storm))$hit)
# Peyd.after.types <- unique(data.frame(tax_table(Peyd.after))$hit)
# 
# Peyd.before.types.subclade <- Peyd.before.types
# Peyd.before.types.subclade <- gsub("_.*","",Peyd.before.types.subclade)
# Peyd.before.types.subclade <- gsub("\\..*","",Peyd.before.types.subclade)
# Peyd.before.types.subclade <- unique(Peyd.before.types.subclade)
# 
# Peyd.storm.types.subclade <- Peyd.storm.types
# Peyd.storm.types.subclade <- gsub("_.*","",Peyd.storm.types.subclade)
# Peyd.storm.types.subclade <- gsub("\\..*","",Peyd.storm.types.subclade)
# Peyd.storm.types.subclade <- unique(Peyd.storm.types.subclade)
# 
# Peyd.after.types.subclade <- Peyd.after.types
# Peyd.after.types.subclade <- gsub("_.*","",Peyd.after.types.subclade)
# Peyd.after.types.subclade <- gsub("\\..*","",Peyd.after.types.subclade)
# Peyd.after.types.subclade <- unique(Peyd.after.types.subclade)

sed.bs <- intersect(sediment.storm.types.subclade,sediment.before.types.subclade)
sed.bsa <- intersect(sed.bs, sediment.after.types.subclade)
sed.bs.only <- length(sed.bs)-length(sed.bsa)
sed.bsa.length <- length(sed.bsa)
sed.ba <- intersect(sediment.before.types.subclade,sediment.after.types.subclade)
sed.ba.only <- length(sed.ba)-length(sed.bsa)
sed.sa <- intersect(sediment.storm.types.subclade,sediment.after.types.subclade)
sed.sa.only <- length(sed.sa)-length(sed.bsa)
sed.a <- length(sediment.after.types.subclade) - sed.sa.only - sed.ba.only - sed.bsa.length
sed.b <- length(sediment.before.types.subclade) - sed.bs.only - sed.ba.only - sed.bsa.length
sed.s <- length(sediment.storm.types.subclade) - sed.sa.only - sed.bs.only - sed.bsa.length


# jpeg(filename="figures/sediment_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
VennDiag <- euler(c("Before" = (sed.b), "Post-Storm" = (sed.s), "After" = (sed.a), "Before&Post-Storm" = (sed.bs.only), "Post-Storm&After" = (sed.sa.only), "Before&After" = (sed.ba.only), "Before&Post-Storm&After" = (sed.bsa.length)))
plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()


water.bs <- intersect(water.storm.types.subclade,water.before.types.subclade)
water.bsa <- intersect(water.bs, water.after.types.subclade)
water.bs.only <- length(water.bs)-length(water.bsa)
water.bsa.length <- length(water.bsa)
water.ba <- intersect(water.before.types.subclade,water.after.types.subclade)
water.ba.only <- length(water.ba)-length(water.bsa)
water.sa <- intersect(water.storm.types.subclade,water.after.types.subclade)
water.sa.only <- length(water.sa)-length(water.bsa)
water.a <- length(water.after.types.subclade) - water.sa.only - water.ba.only - water.bsa.length
water.b <- length(water.before.types.subclade) - water.bs.only - water.ba.only - water.bsa.length
water.s <- length(water.storm.types.subclade) - water.sa.only - water.bs.only - water.bsa.length



# jpeg(filename="figures/water_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
VennDiag2 <- euler(c("Before" = (water.b), "Post-Storm" = (water.s), "After" = (water.a), "Before&Post-Storm" = (water.bs.only), "Post-Storm&After" = (water.sa.only), "Before&After" = (water.ba.only), "Before&Post-Storm&After" = (water.bsa.length)))
plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()

coral.bs <- intersect(coral.storm.types.subclade,coral.before.types.subclade)
coral.bsa <- intersect(coral.bs, coral.after.types.subclade)
coral.bs.only <- length(coral.bs)-length(coral.bsa)
coral.bsa.length <- length(coral.bsa)
coral.ba <- intersect(coral.before.types.subclade,coral.after.types.subclade)
coral.ba.only <- length(coral.ba)-length(coral.bsa)
coral.sa <- intersect(coral.storm.types.subclade,coral.after.types.subclade)
coral.sa.only <- length(coral.sa)-length(coral.bsa)
coral.a <- length(coral.after.types.subclade) - coral.sa.only - coral.ba.only - coral.bsa.length
coral.b <- length(coral.before.types.subclade) - coral.bs.only - coral.ba.only - coral.bsa.length
coral.s <- length(coral.storm.types.subclade) - coral.sa.only - coral.bs.only - coral.bsa.length


# jpeg(filename="figures/coral_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
VennDiag3 <- euler(c("Before" = (coral.b), "Post-Storm" = (coral.s), "After" = (coral.a), "Before&Post-Storm" = (coral.bs.only), "Post-Storm&After" = (coral.sa.only), "Before&After" = (coral.ba.only), "Before&Post-Storm&After" = (coral.bsa.length)))
plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()

# Peyd.bs <- intersect(Peyd.storm.types.subclade,Peyd.before.types.subclade)
# Peyd.bsa <- intersect(Peyd.bs, Peyd.after.types.subclade)
# Peyd.bs.only <- length(Peyd.bs)-length(Peyd.bsa)
# Peyd.bsa.length <- length(Peyd.bsa)
# Peyd.ba <- intersect(Peyd.before.types.subclade,Peyd.after.types.subclade)
# Peyd.ba.only <- length(Peyd.ba)-length(Peyd.bsa)
# Peyd.sa <- intersect(Peyd.storm.types.subclade,Peyd.after.types.subclade)
# Peyd.sa.only <- length(Peyd.sa)-length(Peyd.bsa)
# Peyd.a <- length(Peyd.after.types.subclade) - Peyd.sa.only - Peyd.ba.only - Peyd.bsa.length
# Peyd.b <- length(Peyd.before.types.subclade) - Peyd.bs.only - Peyd.ba.only - Peyd.bsa.length
# Peyd.s <- length(Peyd.storm.types.subclade) - Peyd.sa.only - Peyd.bs.only - Peyd.bsa.length

# jpeg(filename="figures/Peyd_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
# VennDiag4 <- euler(c("Before" = (Peyd.b), "Post-Storm" = (Peyd.s), "After" = (Peyd.a), "Before&Post-Storm" = (Peyd.bs.only), "Post-Storm&After" = (Peyd.sa.only), "Before&After" = (Peyd.ba.only), "Before&Post-Storm&After" = (Peyd.bsa.length)))
# plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
#      fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()



vd1 <- plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = F,main="Sediment")
# vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
#             fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","Post-Storm","After"),main="Water",key=list(space="right",boxes=list(col=timecols), lwd=c(1,4,2),text=list(c("Before","Post-Storm","After"))
#             )) # removed: auto.key=T
vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","Post-Storm","After"),main="Water",auto.key = F)
vd3 <- plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = F,main="Coral")

jpeg(filename="figures/venn_by_fieldseason.jpg", 
     width = 12, height = 4, units="in",res = 300)
grid.arrange(vd1,vd2,vd3,ncol=3)
dev.off()

