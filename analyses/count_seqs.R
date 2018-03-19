# KI Compartment Betadispersion

rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

total_seqs <- sum(taxa_sums(phy97.f.c))
total_seqs

coral_seqs <- sum(taxa_sums(phy97.f.c.coral))
sediment_seqs <- sum(taxa_sums(phy97.f.c.sediment))
water_seqs <- sum(taxa_sums(phy97.f.c.water))
coral_seqs+sediment_seqs+water_seqs
