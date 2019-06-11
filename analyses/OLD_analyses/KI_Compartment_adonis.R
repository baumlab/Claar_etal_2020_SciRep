# KI Compartment adonis

# Clear working directory
rm(list=ls())

# Load necessary files
load("data/KI_Compartment_f_coral_grouped.RData")

# Load necessary libraries
library(vegan)

## Sediment - By field season and site
# Extract sample data as metadata
sediment.metadata <- as(sample_data(phyASV.f.c.sediment), "data.frame")
# Run adonis with wunifrac (weighted unifrac)
adonis(phyloseq::distance(phyASV.f.c.sediment, method="wunifrac") ~ field_season + site,
       data = sediment.metadata)

## Water - By field season and site
water.metadata <- as(sample_data(phyASV.f.c.water), "data.frame")
adonis(distance(phyASV.f.c.water, method="wunifrac") ~ field_season + site,
       data = water.metadata)

## Coral - By field season and site
coral.sp.metadata <- as(sample_data(phyASV.f.c.coral), "data.frame")
adonis(distance(phyASV.f.c.coral, method="wunifrac") ~ Coral_Species + field_season + site, data = coral.sp.metadata)

## Sediment - By field season and disturbance level
adonis(distance(phyASV.f.c.sediment, method="wunifrac") ~ field_season + Dist,
       data = sediment.metadata)
## Water - By field season and disturbance level
adonis(distance(phyASV.f.c.water, method="wunifrac") ~ field_season + Dist,
       data = water.metadata)
## Coral - By field season and disturbance level
adonis(distance(phyASV.f.c.coral, method="wunifrac") ~ Coral_Species + field_season + Dist,data = coral.sp.metadata)
