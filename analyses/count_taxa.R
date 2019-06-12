# KI Compartment - Count sequences
library(phyloseq)

# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

phyASV.f.c.coral.MAeq.p <- transform_sample_counts(phyASV.f.c.coral.MAeq, 
                                                   function(x) x/sum(x))
phyASV.f.c.coral.Peyd.p <- transform_sample_counts(phyASV.f.c.coral.Peyd, 
                                                   function(x) x/sum(x))
phyASV.f.c.coral.Plob.p <- transform_sample_counts(phyASV.f.c.coral.Plob, 
                                                   function(x) x/sum(x))


# Count number of ASVs
ntaxa(phyASV.f.c.coral)
ntaxa(phyASV.f.c.sediment)
ntaxa(phyASV.f.c.water)
ntaxa(phyASV.f.c.coral.MAeq)
ntaxa(phyASV.f.c.coral.Peyd)
ntaxa(phyASV.f.c.coral.Plob)
ntaxa(phyASV.f.c)

# Count number of clades
unique(data.frame(tax_table(phyASV.f.c.coral))$Genus)
unique(data.frame(tax_table(phyASV.f.c.water))$Genus)
unique(data.frame(tax_table(phyASV.f.c.sediment))$Genus)
unique(data.frame(tax_table(phyASV.f.c.coral.MAeq))$Genus)
unique(data.frame(tax_table(phyASV.f.c.coral.Peyd))$Genus)
unique(data.frame(tax_table(phyASV.f.c.coral.Plob))$Genus)
unique(data.frame(tax_table(phyASV.f.c))$Genus)

colMax <- function(X) apply(X, 2, max)
rowMax <- function(X) apply(X, 1, max)

water.ASVs.30plus <- data.frame(otu_table(phyASV.f.c.water.p))[colMax(data.frame(otu_table(phyASV.f.c.water.p)))>0.2499]
water.ASVs.30plus[water.ASVs.30plus<0.24999] <- NA
water.ASVs.dom <- data.frame(colSums(!is.na(water.ASVs.30plus)))
water.ASVs.dom.n <- nrow(water.ASVs.dom)
water.ASVs.dom.names <- rownames(water.ASVs.dom)

phyASV.f.c.sediment.p <- subset_samples(phyASV.f.c.sediment.p, SampleID!="KI15aSSYM029")
sediment.ASVs.30plus <- data.frame(otu_table(phyASV.f.c.sediment.p))[colMax(data.frame(otu_table(phyASV.f.c.sediment.p)))>0.2499]
sediment.ASVs.30plus[sediment.ASVs.30plus<0.24999] <- NA
sediment.ASVs.dom <- data.frame(colSums(!is.na(sediment.ASVs.30plus)))
sediment.ASVs.dom.n <- nrow(sediment.ASVs.dom)
sediment.ASVs.dom.names <- rownames(sediment.ASVs.dom)

MAeq.ASVs.30plus <- data.frame(otu_table(phyASV.f.c.coral.MAeq.p))[colMax(data.frame(otu_table(phyASV.f.c.coral.MAeq.p)))>0.2499]
MAeq.ASVs.30plus[MAeq.ASVs.30plus<0.24999] <- NA
MAeq.ASVs.dom <- data.frame(colSums(!is.na(MAeq.ASVs.30plus)))
MAeq.ASVs.dom.n <- nrow(MAeq.ASVs.dom)
MAeq.ASVs.dom.names <- rownames(MAeq.ASVs.dom)

Peyd.ASVs.30plus <- data.frame(otu_table(phyASV.f.c.coral.Peyd.p))[colMax(data.frame(otu_table(phyASV.f.c.coral.Peyd.p)))>0.2499]
Peyd.ASVs.30plus[Peyd.ASVs.30plus<0.24999] <- NA
Peyd.ASVs.dom <- data.frame(colSums(!is.na(Peyd.ASVs.30plus)))
Peyd.ASVs.dom.n <- nrow(Peyd.ASVs.dom)
Peyd.ASVs.dom.names <- rownames(Peyd.ASVs.dom)

Plob.ASVs.30plus <- data.frame(otu_table(phyASV.f.c.coral.Plob.p))[colMax(data.frame(otu_table(phyASV.f.c.coral.Plob.p)))>0.2499]
Plob.ASVs.30plus[Plob.ASVs.30plus<0.24999] <- NA
Plob.ASVs.dom <- data.frame(colSums(!is.na(Plob.ASVs.30plus)))
Plob.ASVs.dom.n <- nrow(Plob.ASVs.dom)
Plob.ASVs.dom.names <- rownames(Plob.ASVs.dom)

water.ASVs.dom.n
sediment.ASVs.dom.n
MAeq.ASVs.dom.n
Peyd.ASVs.dom.n
Plob.ASVs.dom.n

water.ASVs.dom.names
sediment.ASVs.dom.names
MAeq.ASVs.dom.names
Peyd.ASVs.dom.names
Plob.ASVs.dom.names

writeXStringSet(file="analyses/refseqs/Plob_dom.fasta",refseq(phyASV.f.c)[Plob.ASVs.dom.names])
writeXStringSet(file="analyses/refseqs/Maeq_dom.fasta",refseq(phyASV.f.c)[MAeq.ASVs.dom.names])
writeXStringSet(file="analyses/refseqs/Peyd_dom.fasta",refseq(phyASV.f.c)[Peyd.ASVs.dom.names])




writeXStringSet(file="analyses/refseqs/ASV3.fasta",refseq(phyASV.f.c)["ASV3"])
writeXStringSet(file="analyses/refseqs/ASV7.fasta",refseq(phyASV.f.c)["ASV7"])
writeXStringSet(file="analyses/refseqs/ASV16.fasta",refseq(phyASV.f.c)["ASV16"])
writeXStringSet(file="analyses/refseqs/ASV28.fasta",refseq(phyASV.f.c)["ASV28"])
writeXStringSet(file="analyses/refseqs/ASV36.fasta",refseq(phyASV.f.c)["ASV36"])
writeXStringSet(file="analyses/refseqs/ASV51.fasta",refseq(phyASV.f.c)["ASV51"])
writeXStringSet(file="analyses/refseqs/ASV105.fasta",refseq(phyASV.f.c)["ASV105"])
