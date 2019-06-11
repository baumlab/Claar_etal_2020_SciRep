# KI Compartment - Count sequences

# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Calculate the total number of sequences
total_seqs <- sum(taxa_sums(phyASV.f.c))
total_seqs

# Calculate the number of sequences in each compartment
coral_seqs <- sum(taxa_sums(phyASV.f.c.coral))
sediment_seqs <- sum(taxa_sums(phyASV.f.c.sediment))
water_seqs <- sum(taxa_sums(phyASV.f.c.water))

# Check that total_seqs equals the sum of all compartments combined
coral_seqs+sediment_seqs+water_seqs
