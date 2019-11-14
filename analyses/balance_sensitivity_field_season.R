# Clear working environment
rm(list=ls())

# Load necessary packages
library(phyloseq)
library(dplyr)
library(vegan)
library(purrr)

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")

# I'm sorry this is very slow. And also that it's not a function.

physeq<-phyASV.f.c.sediment
nreps<-100
nbal<-22

# subset_balanced <- function(physeq,nreps,nbal) {
ado.F.list <- list()
ado.R2.list <- list()
ado.P.list <- list()

for(n in 1:nreps){
  subbal <- sample_data(physeq) %>% 
    group_by(field_season) %>% 
    sample_n(size = nbal)
  # nrow(subbal[subbal$Dist=="VeryHigh",])
  # nrow(subbal[subbal$Dist=="Medium",])
  
  phy_sub <- subset_samples(physeq,
                            sample_data(physeq)$SampleID %in% subbal$SampleID)
  # nrow(sample_data(phy_sub)[sample_data(phy_sub)$Dist=="VeryHigh",])
  # nrow(sample_data(phy_sub)[sample_data(phy_sub)$Dist=="Medium",])
  
  set.seed(2020)
  ufdist <- UniFrac(phy_sub, weighted=T, 
                    normalized=F, parallel=F, fast=T)
  # sediment.bd.dist <- betadisper(d=ufdist, 
  #                                group=sample_data(phy_sub)$Dist,
  #                                type="centroid", bias.adjust=FALSE)
  metadata <- as(sample_data(phy_sub), "data.frame")
  ado <- adonis(ufdist ~ Dist, data = metadata)
  set.seed(NULL)
  
  ado.F.list[[n]] <- ado$aov.tab$F.Model[1]
  ado.R2.list[[n]] <- ado$aov.tab$R2[1]
  ado.P.list[[n]] <- ado$aov.tab$`Pr(>F)`[1]
  print(n)
}
ado.F <- do.call(rbind,ado.F.list)
ado.R2 <- do.call(rbind,ado.R2.list)
ado.P <- do.call(rbind,ado.P.list)
ado.stats <- do.call(cbind, list(ado.F, ado.R2, ado.P))
# }
# 
# subset_balanced(phyASV.f.c.sediment,3,33)
ado.stats
