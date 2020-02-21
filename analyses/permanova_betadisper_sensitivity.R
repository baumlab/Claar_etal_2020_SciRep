# KI Compartment Betadispersion - Sensitivity

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")

library(vegan)
library(ggplot2)
library(phyloseq)

phyASV.f.c.coral.Peyd <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phyASV.f.c.coral.MAeq <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Montipora_foliosa")
phyASV.f.c.coral.Plob <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Porites_lobata")

phyASV.f.c.M <- subset_samples(phyASV.f.c,sample_data(phyASV.f.c)$Dist=="Medium")
phyASV.f.c.VH <- subset_samples(phyASV.f.c,sample_data(phyASV.f.c)$Dist=="VeryHigh")

phyASV.f.c.sediment.2014 <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$field_season=="KI2014")
phyASV.f.c.sediment.2015a <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$field_season=="KI2015a_Post")
phyASV.f.c.sediment.2015b <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$field_season=="KI2015b")


nrow(sample_data(phyASV.f.c.sediment.M)) #33
nrow(sample_data(phyASV.f.c.sediment.VH)) #34
nrow(sample_data(phyASV.f.c.water.M)) #25
nrow(sample_data(phyASV.f.c.water.VH)) #25
# nrow(sample_data(phyASV.f.c.coral.MAeq.M)) #40 # dont need adnois is insig.
# nrow(sample_data(phyASV.f.c.coral.MAeq.VH)) #38
nrow(sample_data(phyASV.f.c.coral.Peyd.M)) #39
nrow(sample_data(phyASV.f.c.coral.Peyd.VH)) #31
# nrow(sample_data(phyASV.f.c.coral.Plob.M)) #34 # dont need, both are insig.
# nrow(sample_data(phyASV.f.c.coral.Plob.VH)) #39
nrow(sample_data(phyASV.f.c.M)) #171
nrow(sample_data(phyASV.f.c.VH)) #167

nrow(sample_data(phyASV.f.c.sediment.2014))
nrow(sample_data(phyASV.f.c.sediment.2015a))
nrow(sample_data(phyASV.f.c.sediment.2015b))

