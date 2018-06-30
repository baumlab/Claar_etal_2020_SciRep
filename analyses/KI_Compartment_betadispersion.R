# KI Compartment Prepare for Betadispersion Analyses and Figures

# Clear working environment
rm(list=ls())

# Load in necessary RData files
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Load necessary packages
library(vegan)

# Subset coral species
phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.Peyd <- prune_taxa(taxa_sums(phy97.f.c.coral.Peyd)>0,phy97.f.c.coral.Peyd)
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.MAeq <- prune_taxa(taxa_sums(phy97.f.c.coral.MAeq)>0,phy97.f.c.coral.MAeq)
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")
phy97.f.c.coral.Plob <- prune_taxa(taxa_sums(phy97.f.c.coral.Plob)>0,phy97.f.c.coral.Plob)

# Calculate Unifrac distances for each compartment
sediment.ufdist <- UniFrac(phy97.f.c.sediment, weighted=T, 
                           normalized=F, parallel=F, fast=T)
water.ufdist <- UniFrac(phy97.f.c.water, weighted=T, 
                        normalized=F, parallel=F, fast=T)
all.ufdist <- UniFrac(phy97.f.c, weighted=T, 
                      normalized=F, parallel=F, fast=T)
Peyd.ufdist <- UniFrac(phy97.f.c.coral.Peyd, weighted=T, 
                       normalized=F, parallel=F, fast=T)
MAeq.ufdist <- UniFrac(phy97.f.c.coral.MAeq, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Plob.ufdist <- UniFrac(phy97.f.c.coral.Plob, weighted=T, 
                       normalized=F, parallel=F, fast=T)
coral.ufdist <- UniFrac(phy97.f.c.coral, weighted=T, 
                        normalized=F, parallel=F, fast=T)

# Run betadisper for each compartment
sediment.bd.dist <- betadisper(d=sediment.ufdist, 
                               group=sample_data(phy97.f.c.sediment)$Dist,
                               type="centroid", bias.adjust=FALSE)
water.bd.dist <- betadisper(d=water.ufdist, 
                            group=sample_data(phy97.f.c.water)$Dist,
                            type="centroid", bias.adjust=FALSE)
all.bd.dist <- betadisper(d=all.ufdist, 
                          group=sample_data(phy97.f.c)$Dist,
                          type="centroid", bias.adjust=FALSE)
Peyd.bd.dist <- betadisper(d=Peyd.ufdist, 
                           group=sample_data(phy97.f.c.coral.Peyd)$Dist,
                           type="centroid", bias.adjust=FALSE)
MAeq.bd.dist <- betadisper(d=MAeq.ufdist, 
                           group=sample_data(phy97.f.c.coral.MAeq)$Dist,
                           type="centroid", bias.adjust=FALSE)
Plob.bd.dist <- betadisper(d=Plob.ufdist, 
                           group=sample_data(phy97.f.c.coral.Plob)$Dist,
                           type="centroid", bias.adjust=FALSE)
coral.bd.dist <- betadisper(d=coral.ufdist, 
                            group=sample_data(phy97.f.c.coral)$Dist,
                            type="centroid", bias.adjust=FALSE)
coral.bd.species <- betadisper(d=coral.ufdist, 
                               group=sample_data(phy97.f.c.coral)$Coral_Species,
                               type="centroid", bias.adjust=FALSE)

# Test betadisper results
betadisper.sediment <- anova(sediment.bd.dist)
betadisper.water <- anova(water.bd.dist)
betadisper.Peyd <- anova(Peyd.bd.dist)
betadisper.MAeq <- anova(MAeq.bd.dist)
betadisper.Plob <- anova(Plob.bd.dist)
betadisper.all <- anova(all.bd.dist)

# Save 
save.image(file="analyses/betadisper.RData")