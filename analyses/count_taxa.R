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

sum(water.ASVs %in% sediment.ASVs)
sum(sediment.ASVs %in% water.ASVs)
sum(sediment.ASVs %in% coral.ASVs)
sum(coral.ASVs %in% sediment.ASVs)
sum(water.ASVs %in% coral.ASVs)
sum(coral.ASVs %in% water.ASVs)

nrow(data.frame(intersect(intersect(water.ASVs,sediment.ASVs),coral.ASVs)))
nrow(data.frame(intersect(water.ASVs,sediment.ASVs)))
nrow(data.frame(intersect(water.ASVs,coral.ASVs)))
nrow(data.frame(intersect(coral.ASVs,sediment.ASVs)))
nrow(data.frame(intersect(intersect(water.ASVs,sediment.ASVs),coral.ASVs)))

outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

water.sed <- c(water.ASVs,sediment.ASVs)
water.coral <- c(water.ASVs,coral.ASVs)
sed.coral <- c(sediment.ASVs,coral.ASVs)
all.ASVs

nrow(data.frame(outersect(water.sed,all.ASVs)))
nrow(data.frame(outersect(water.coral,all.ASVs)))
nrow(data.frame(outersect(sed.coral,all.ASVs)))
