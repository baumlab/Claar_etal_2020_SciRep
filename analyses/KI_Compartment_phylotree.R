rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")

phy97.f.c.sediment <- subset_samples(phy97.f.c.sediment,field_season!="KI2015c")

col <- c(KI2014 = "#2c7fb8", KI2015a_Pre = "#7fcdbb", KI2015a_Post = "#253494", KI2015b = "#41b6c4")
cols2 <- c(water = "#5ab4ac", sediment = "#d8b365")

water <- subset_taxa(physeq = phy97.f.c.water,taxa_sums(phy97.f.c.water) > 10)
sediment <- subset_taxa(physeq = phy97.f.c.sediment,taxa_sums(phy97.f.c.sediment) > 10)

pt1 <- plot_tree(water,color="field_season",label.tips = "hit") + 
  scale_color_manual(values=col) + 
  guides(fill=F,color=F) + 
  ggtitle("Water")
pt2 <- plot_tree(sediment,color="field_season", label.tips = "hit") + 
  scale_color_manual(values=col) + 
  ggtitle("Sediment")

jpeg(filename="figures/sediment_water_phylotrees.jpg", width = 11, height = 8, units="in",res = 300)
grid.arrange(pt1,pt2,ncol=2)
dev.off()
