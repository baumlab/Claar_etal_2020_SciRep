# KI Compartment adonis

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
library(vegan)

sediment.metadata <- as(sample_data(phy97.f.c.sediment), "data.frame")
adonis(distance(phy97.f.c.sediment, method="wunifrac") ~ field_season + site,
       data = sediment.metadata)

water.metadata <- as(sample_data(phy97.f.c.water), "data.frame")
adonis(distance(phy97.f.c.water, method="wunifrac") ~ field_season + site,
       data = water.metadata)

coral.sp.metadata <- as(sample_data(phy97.f.c.coral), "data.frame")
adonis(distance(phy97.f.c.coral, method="wunifrac") ~ Coral_Species + field_season + site,
       data = coral.sp.metadata)

adonis(distance(phy97.f.c.sediment, method="wunifrac") ~ field_season + Dist,
       data = sediment.metadata)
adonis(distance(phy97.f.c.water, method="wunifrac") ~ field_season + Dist,
       data = water.metadata)
adonis(distance(phy97.f.c.coral, method="wunifrac") ~ Coral_Species + field_season + Dist,
       data = coral.sp.metadata)
