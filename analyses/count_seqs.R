# KI Compartment - Count sequences

# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Calculate the total number of sequences
total_seqs <- sum(taxa_sums(phy97.f.c))
total_seqs

# Calculate the number of sequences in each compartment
coral_seqs <- sum(taxa_sums(phy97.f.c.coral))
sediment_seqs <- sum(taxa_sums(phy97.f.c.sediment))
water_seqs <- sum(taxa_sums(phy97.f.c.water))

# Check that total_seqs equals the sum of all compartments combined
coral_seqs+sediment_seqs+water_seqs
