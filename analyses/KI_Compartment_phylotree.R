rm(list=ls())

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")


phy97.f.c.sediment <- subset_samples(phy97.f.c.sediment,field_season!="KI2015c")

water <- subset_taxa(physeq = phy97.f.c.water,taxa_sums(phy97.f.c.water) > 10)
sediment <- subset_taxa(physeq = phy97.f.c.sediment,taxa_sums(phy97.f.c.sediment) > 10)

pt1 <- plot_tree(water,color="field_season",label.tips = "hit") + 
  scale_color_manual(values=timecols) + 
  guides(fill=F,color=F) + 
  ggtitle("Water")
pt2 <- plot_tree(sediment,color="field_season", label.tips = "hit") + 
  scale_color_manual(values=timecols) + 
  ggtitle("Sediment")

jpeg(filename="figures/sediment_water_phylotrees.jpg", width = 11, height = 8, units="in",res = 300)
grid.arrange(pt1,pt2,ncol=2)
dev.off()
