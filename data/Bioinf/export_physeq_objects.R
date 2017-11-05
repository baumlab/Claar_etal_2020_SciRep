# Exporting phyloseq objects

load("analyses/KI_Compartment.RData")

# Write tax table for phylogenetically-informed diversity
write.table(data.frame(phyloseq::tax_table(phy.f)), "data/tax_table.txt", row.names=T, quote=F)
# Write otu table for phylogenetically-informed diversity
write.table(data.frame(phyloseq::otu_table(phy.f)), "data/otu_table.txt", row.names=T, quote=F)
