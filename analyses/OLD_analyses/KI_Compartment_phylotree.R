# Load necessary packages
library(ggplot2)

# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Subset taxa, remove any with less than 10 reads
water <- subset_taxa(physeq = phyASV.f.c.water,taxa_sums(phyASV.f.c.water) > 10)
sediment <- subset_taxa(physeq = phyASV.f.c.sediment,taxa_sums(phyASV.f.c.sediment) > 10)

# Plot tree - Water
pt1 <- plot_tree(water,color="field_season",label.tips = "Genus") + 
  scale_color_manual(values=timecols) + 
  guides(fill=F,color=F) + 
  ggtitle("Water")
# Plot tree - Sediment
pt2 <- plot_tree(sediment,color="field_season", label.tips = "Genus") + 
  scale_color_manual(values=timecols) + 
  ggtitle("Sediment")

# Create jpg
jpeg(filename="figures/sediment_water_phylotrees.jpg", width = 11, height = 8, units="in",res = 300)
grid.arrange(pt1,pt2,ncol=2)
dev.off()
