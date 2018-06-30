# Clear the working environment
rm(list=ls())

# Load necessary packages
library(vegan)

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Helper function from http://joey711.github.io/phyloseq-demo/phyloseq-demo.html to convert phyloseq objects to vegan-friendly formatting
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Subset taxa, removing any taxa that have fewer than 50 reads
water <- subset_taxa(physeq = phy97.f.c.water,taxa_sums(phy97.f.c.water) > 50)
sediment <- subset_taxa(physeq = phy97.f.c.sediment,taxa_sums(phy97.f.c.sediment) > 50)
all <- subset_taxa(physeq = phy97.f.c,taxa_sums(phy97.f.c) > 50)

# Transform to vegan-friendly objects
water.veg <- veganotu(water)
sediment.veg <- veganotu(sediment)
all.veg <- veganotu(all)

# Indicator species analysis - Water by Field Season
groups.water <- as.numeric(sample_data(water)$field_season)
indval.water = multipatt(as.data.frame(water.veg), groups.water, control = how(nperm=999))
summary(indval.water)

# Indicator species analysis - Sediment by Field Season
groups.sed <- as.numeric(sample_data(sediment)$field_season)
indval.sed = multipatt(as.data.frame(sediment.veg), groups.sed, control = how(nperm=999))
summary(indval.sed)

# Indicator species analysis - All by Field Season
groups.all <- as.numeric(sample_data(all)$SampleType)
indval.all = multipatt(as.data.frame(all.veg), groups.all, control = how(nperm=999))
summary(indval.all)
