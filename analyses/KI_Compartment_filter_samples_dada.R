# Import Libraries
library(stringr)
library(reshape2)
library(phyloseq)
library(seqinr)
library(phangorn) 
library(caroline)
library(DECIPHER)


# Clear workspace
rm(list=ls())

# Load filtered RData object from output of filter_notsym.R script
load("analyses/KI_Compartment_dada.RData")

sample_data(phy.f)$field_season <- sample_data(phy.f)$Year_Pre_Post
sample_data(phy.f)$field_season <- as.character(sample_data(phy.f)$field_season)
sample_data(phy.f)$field_season <- replace(sample_data(phy.f)$field_season, sample_data(phy.f)$field_season=="2014", "KI2014")
sample_data(phy.f)$field_season <- replace(sample_data(phy.f)$field_season, sample_data(phy.f)$field_season=="2015Jan_Pre", "KI2015a_Pre")
sample_data(phy.f)$field_season <- replace(sample_data(phy.f)$field_season, sample_data(phy.f)$field_season=="2015Jan_Post", "KI2015a_Post")
sample_data(phy.f)$field_season <- replace(sample_data(phy.f)$field_season, sample_data(phy.f)$field_season=="2015May", "KI2015b")
sample_data(phy.f)$field_season <- replace(sample_data(phy.f)$field_season, sample_data(phy.f)$field_season=="2015July", "KI2015c")

colnames(sample_data(phy.f))[colnames(sample_data(phy.f))=="Site"] <- "site"

################################### Filtering ######################################
# Remove samples that were sequenced with this set, but are not included in this ms
# Remove KI2015c samples
phy.f <- prune_samples(sample_data(phy.f)$field_season!="KI2015c", phy.f)
# Remove Site34 samples
phy.f <- prune_samples(sample_data(phy.f)$site!="34", phy.f)
phy.f <- prune_taxa(taxa_sums(phy.f)>0, phy.f)

############################## Site Formatting ####################################

# Characterize sites by disturbance level
VeryHigh <- c(30,31,32)
High <- c(1,6,25,26,38,40)
Medium <- c(7,8,12,13,14,22,33,34,35)
Low <- c(2,3,4,9,23,24)
VeryLow <- c(5,10,11,15,16,17,18,19,20,21,28,29,36,37,39)

# Make site a factor
sample_data(phy.f)$site <- factor(sample_data(phy.f)$site)
# Create a disturbance level column, start as site
sample_data(phy.f)$Dist <- sample_data(phy.f)$site

# Finish creating Dist
for (i in VeryHigh){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryHigh"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("VeryHigh")))
}

for (i in High){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("High"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("High")))
}

for (i in Medium){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Medium"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("Medium")))
}

for (i in Low){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Low"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("Low")))
}

for (i in VeryLow){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryLow"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("VeryLow")))
}

# Order levels
levels(sample_data(phy.f)$Dist) <- c("VeryHigh","High","Medium","Low","VeryLow")

###################### Physeq formatting and tree #####################################
# Assign new name for clarity
phyASV.f.c <- phy.f
phyASV.f.c <- prune_taxa(taxa_sums(phyASV.f.c)>0,phyASV.f.c)

# Write tax table for phylogenetically-informed diversity
write.table(data.frame(tax_table(phyASV.f.c)), "data/tax_table.txt", row.names=T, quote=F)
# Write otu table for phylogenetically-informed diversity
write.delim(data.frame(otu_table(phyASV.f.c)), "data/otu_table.tsv", quote = FALSE, row.names = T, sep = "\t")

ASV.As <- subset_taxa(phyASV.f.c,Genus=="g__Symbiodinium")
ASV.Bs <- subset_taxa(phyASV.f.c,Genus=="g__Breviolum")
ASV.Cs <- subset_taxa(phyASV.f.c,Genus=="g__Cladocopium")
ASV.Ds <- subset_taxa(phyASV.f.c,Genus=="g__Durusdinium")
ASV.Es <- subset_taxa(phyASV.f.c,Genus=="g__Effrenium")
ASV.Fs <- subset_taxa(phyASV.f.c,Genus=="g__Fugacium")
ASV.Gs <- subset_taxa(phyASV.f.c,Genus=="g__Gerakladium")
# ASV.Hs <- subset_taxa(phyASV.f.c,Genus=="g__cladeH") # no seqs
ASV.Is <- subset_taxa(phyASV.f.c,Genus=="g__cladeI")

ASV.seqs <- refseq(phyASV.f.c)
min(width(ASV.seqs))
max(width(ASV.seqs))

ASV.A.seqs <- refseq(ASV.As)
ASV.A.seqs.aligned <- AlignSeqs(ASV.A.seqs)
writeXStringSet(ASV.A.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_A_tree_seqs_aligned.fasta")
ASV.B.seqs <- refseq(ASV.Bs)
ASV.B.seqs.aligned <- refseq(ASV.Bs) # There is only one sequence here! Can't align
writeXStringSet(ASV.B.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_B_tree_seqs_aligned.fasta")
ASV.C.seqs <- refseq(ASV.Cs)
ASV.C.seqs.aligned <- AlignSeqs(ASV.C.seqs)
writeXStringSet(ASV.C.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_C_tree_seqs_aligned.fasta")
ASV.D.seqs <- refseq(ASV.Ds)
ASV.D.seqs.aligned <- AlignSeqs(ASV.D.seqs)
writeXStringSet(ASV.D.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_D_tree_seqs_aligned.fasta")
ASV.E.seqs <- refseq(ASV.Es)
writeXStringSet(ASV.E.seqs, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_Es.fasta")
ASV.E.seqs.aligned <- AlignSeqs(ASV.E.seqs)
writeXStringSet(ASV.E.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_E_tree_seqs_aligned.fasta")
ASV.F.seqs <- refseq(ASV.Fs)
ASV.F.seqs.aligned <- AlignSeqs(ASV.F.seqs)
writeXStringSet(ASV.F.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_F_tree_seqs_aligned.fasta")
ASV.G.seqs <- refseq(ASV.Gs)
ASV.G.seqs.aligned <- AlignSeqs(ASV.G.seqs)
writeXStringSet(ASV.G.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_G_tree_seqs_aligned.fasta")
ASV.I.seqs <- refseq(ASV.Is)
ASV.I.seqs.aligned <- AlignSeqs(ASV.I.seqs)
writeXStringSet(ASV.I.seqs.aligned, # write the alignment to a new FASTA file
                file="data/Bioinf/tree/ASV_I_tree_seqs_aligned.fasta")

min(width(ASV.A.seqs))
max(width(ASV.A.seqs))
min(width(ASV.B.seqs))
max(width(ASV.B.seqs))
min(width(ASV.C.seqs))
max(width(ASV.C.seqs))
min(width(ASV.D.seqs))
max(width(ASV.D.seqs))
min(width(ASV.E.seqs))
max(width(ASV.E.seqs))
min(width(ASV.F.seqs))
max(width(ASV.F.seqs))
min(width(ASV.G.seqs))
max(width(ASV.G.seqs))
min(width(ASV.I.seqs))
max(width(ASV.I.seqs))



#https://rdrr.io/rforge/seqinr/man/dist.alignment.html
#returns sqrt of pairwise genetic distance, then squared the matrices
A.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_A_tree_seqs_aligned.fasta", format= "fasta")
A.dis <- (as.matrix(dist.alignment(A.seqs, matrix = "identity" )))^2
write.csv(A.dis, file="data/Bioinf/tree/ASV_A_dis_matx.csv")

B.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_B_tree_seqs_aligned.fasta", format= "fasta")
B.dis <- (as.matrix(dist.alignment(B.seqs, matrix = "identity" )))^2
write.csv(B.dis, file="data/Bioinf/tree/ASV_B_dis_matx.csv")

C.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_C_tree_seqs_aligned.fasta", format= "fasta")
C.dis <- (as.matrix(dist.alignment(C.seqs, matrix = "identity" )))^2
write.csv(C.dis, file="data/Bioinf/tree/ASV_C_dis_matx.csv")

D.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_D_tree_seqs_aligned.fasta", format= "fasta")
D.dis <- (as.matrix(dist.alignment(D.seqs, matrix = "identity" )))^2
write.csv(D.dis, file="data/Bioinf/tree/ASV_D_dis_matx.csv")

E.seqs <- read.alignment(file = "data/BioinE/tree/ASV_E_tree_seqs_aligned.fasta", format= "fasta")
E.dis <- (as.matrix(dist.alignment(E.seqs, matrix = "identity" )))^2
write.csv(E.dis, file="data/Bioinf/tree/ASV_E_dis_matx.csv")

F.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_F_tree_seqs_aligned.fasta", format= "fasta")
F.dis <- (as.matrix(dist.alignment(F.seqs, matrix = "identity" )))^2
write.csv(F.dis, file="data/Bioinf/tree/ASV_F_dis_matx.csv")

G.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_G_tree_seqs_aligned.fasta", format= "fasta")
G.dis <- (as.matrix(dist.alignment(G.seqs, matrix = "identity" )))^2
write.csv(G.dis, file="data/Bioinf/tree/ASV_G_dis_matx.csv")

I.seqs <- read.alignment(file = "data/Bioinf/tree/ASV_I_tree_seqs_aligned.fasta", format= "fasta")
I.dis <- (as.matrix(dist.alignment(I.seqs, matrix = "identity" )))^2
write.csv(I.dis, file="data/Bioinf/tree/ASV_I_dis_matx.csv")

#give clade distances using average 28s distance from Pochon and Gates 2010
A_B <- matrix(0.219, ncol=ncol(A.dis), nrow=nrow(B.dis), dimnames=list(rownames(B.dis), colnames(A.dis)))
A_C <- matrix(0.1960, ncol=ncol(A.dis), nrow=nrow(C.dis), dimnames=list(rownames(C.dis), colnames(A.dis)))
A_D <- matrix(0.1775, ncol=ncol(A.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(A.dis)))
A_E <- matrix(0.2085, ncol=ncol(A.dis), nrow=nrow(E.dis), dimnames=list(rownames(E.dis), colnames(A.dis)))
A_F <- matrix(0.2085, ncol=ncol(A.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(A.dis)))
A_G <- matrix(0.216, ncol=ncol(A.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(A.dis)))
A_I <- matrix(0.205, ncol=ncol(A.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(A.dis)))
B_C <- matrix(0.114, ncol=ncol(B.dis), nrow=nrow(C.dis), dimnames=list(rownames(C.dis), colnames(B.dis)))
B_D <- matrix(0.1705, ncol=ncol(B.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(B.dis)))
B_E <- matrix(0.1355, ncol=ncol(B.dis), nrow=nrow(E.dis), dimnames=list(rownames(E.dis), colnames(B.dis)))
B_F <- matrix(0.1355, ncol=ncol(B.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(B.dis)))
B_G <- matrix(0.21, ncol=ncol(B.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(B.dis)))
B_I <- matrix(0.151, ncol=ncol(B.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(B.dis)))
C_D <- matrix(0.1520, ncol=ncol(C.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(C.dis)))
C_E <- matrix(0.0945, ncol=ncol(C.dis), nrow=nrow(E.dis), dimnames=list(rownames(E.dis), colnames(C.dis)))
C_F <- matrix(0.0945, ncol=ncol(C.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(C.dis)))
C_G <- matrix(0.187, ncol=ncol(C.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(C.dis)))
C_I <- matrix(0.137, ncol=ncol(C.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(C.dis)))
D_E <- matrix(0.1691, ncol=ncol(D.dis), nrow=nrow(E.dis), dimnames=list(rownames(E.dis), colnames(D.dis)))
D_F <- matrix(0.1691, ncol=ncol(D.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(D.dis)))
D_G <- matrix(0.1795, ncol=ncol(D.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(D.dis)))
D_I <- matrix(0.169125, ncol=ncol(D.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(D.dis)))
E_F <- matrix(0.1691, ncol=ncol(E.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(E.dis)))
E_G <- matrix(0.1691, ncol=ncol(E.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(E.dis)))
E_I <- matrix(0.1691, ncol=ncol(E.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(E.dis)))
F_G <- matrix(0.2072, ncol=ncol(F.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(F.dis)))
F_I <- matrix(0.14925, ncol=ncol(F.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(F.dis)))
G_I <- matrix(0.194, ncol=ncol(G.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(G.dis)))

#build ACDG matrix
col1 <- rbind(A.dis, A_B, A_C, A_D, A_F, A_G, A_I)
col2 <- rbind(matrix(NA, nrow=nrow(A.dis), ncol=ncol(B.dis), dimnames=list(rownames(A.dis), colnames(B.dis))), B.dis, B_C, B_D, B_F, B_G, B_I)
col3 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis), ncol=ncol(C.dis), dimnames=list(c(rownames(A.dis), rownames(B.dis)), colnames(C.dis))), C.dis, C_D, C_F, C_G, C_I)
col4 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis), ncol=ncol(D.dis), dimnames=list(c(rownames(A.dis), rownames(B.dis), rownames(C.dis)), colnames(D.dis))), D.dis, D_F, D_G, D_I)
col5 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis)+nrow(D.dis), ncol=ncol(F.dis), dimnames=list(c(rownames(A.dis), rownames(B.dis), rownames(C.dis), rownames(D.dis)), colnames(F.dis))), F.dis, F_G, F_I)
col6 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis)+nrow(D.dis)+nrow(F.dis), ncol=ncol(G.dis), dimnames=list(c(rownames(A.dis), rownames(B.dis), rownames(C.dis), rownames(D.dis),  rownames(F.dis)), colnames(G.dis))), G.dis, G_I)
col7 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis)+nrow(D.dis)+nrow(F.dis)+nrow(G.dis), ncol=ncol(I.dis), dimnames=list(c(rownames(A.dis), rownames(B.dis), rownames(C.dis), rownames(D.dis),  rownames(F.dis), rownames(G.dis)), colnames(I.dis))), I.dis)

ubermatrix <- cbind(col1, col2, col3, col4, col5, col6, col7)
dim(ubermatrix)

#build tree
uber.tree <- phangorn::upgma(ubermatrix)
plot(uber.tree, main="UPGMA")

#write tree to file
write.tree(uber.tree, file="data/Bioinf/tree/uber.tre")

#use tree and OTU table for calculating beta_diversity.py
#http://qiime.org/scripts/beta_diversity.html

class(uber.tree)

# Slot uber tree into the phy_tree slot of the phyloseq object
phy_tree(phyASV.f.c) <- phy_tree(uber.tree)

######################### Subset by compartment and calc seqs ##################################

# Transform sample counts to proportional abundance for downstream analyses
phyASV.f.c.p <- transform_sample_counts(phyASV.f.c, function(x) x/sum(x))

# Note, this is now the same as phyASV.f.c - renaming and subsetting was completed before when water and sediment samples were included
phyASV.f.c.coral <- subset_samples(phyASV.f.c,SampleType=="coral")
phyASV.f.c.coral <- subset_taxa(phyASV.f.c.coral, taxa_sums(phyASV.f.c.coral) > 0, prune=TRUE)

phyASV.f.c.water <- subset_samples(phyASV.f.c,SampleType=="water")
phyASV.f.c.water <- subset_taxa(phyASV.f.c.water, taxa_sums(phyASV.f.c.water) > 0, prune=TRUE)

phyASV.f.c.sediment <- subset_samples(phyASV.f.c,SampleType=="sediment")
phyASV.f.c.sediment <- subset_taxa(phyASV.f.c.sediment, taxa_sums(phyASV.f.c.sediment) > 0, prune=TRUE)

# Transform sample counts to proportional abundance for downstream analyses
phyASV.f.c.coral.p <- transform_sample_counts(phyASV.f.c.coral, function(x) x/sum(x))
phyASV.f.c.water.p <- transform_sample_counts(phyASV.f.c.water, function(x) x/sum(x))
phyASV.f.c.sediment.p <- transform_sample_counts(phyASV.f.c.sediment, function(x) x/sum(x))

# Calculate number of sequences in phyASV.f.c
total_seqs <- sum(taxa_sums(phyASV.f.c))
coral_seqs <- sum(taxa_sums(phyASV.f.c.coral))
water_seqs <- sum(taxa_sums(phyASV.f.c.water))
sediment_seqs <- sum(taxa_sums(phyASV.f.c.sediment))

############################## Subset by compartment/disturbance #######################
# Subset by compartment and disturbance level
phyASV.f.c.coral.VH <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Dist=="VeryHigh")
phyASV.f.c.coral.VH <- prune_taxa(taxa_sums(phyASV.f.c.coral.VH)>0,phyASV.f.c.coral.VH)
phyASV.f.c.coral.M <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Dist=="Medium")
phyASV.f.c.coral.M <- prune_taxa(taxa_sums(phyASV.f.c.coral.M)>0,phyASV.f.c.coral.M)
phyASV.f.c.sediment.VH <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$Dist=="VeryHigh")
phyASV.f.c.sediment.VH <- prune_taxa(taxa_sums(phyASV.f.c.sediment.VH)>0,phyASV.f.c.sediment.VH)
phyASV.f.c.sediment.M <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$Dist=="Medium")
phyASV.f.c.sediment.M <- prune_taxa(taxa_sums(phyASV.f.c.sediment.M)>0,phyASV.f.c.sediment.M)
phyASV.f.c.water.VH <- subset_samples(phyASV.f.c.water,sample_data(phyASV.f.c.water)$Dist=="VeryHigh")
phyASV.f.c.water.VH <- prune_taxa(taxa_sums(phyASV.f.c.water.VH)>0,phyASV.f.c.water.VH)
phyASV.f.c.water.M <- subset_samples(phyASV.f.c.water,sample_data(phyASV.f.c.water)$Dist=="Medium")
phyASV.f.c.water.M <- prune_taxa(taxa_sums(phyASV.f.c.water.M)>0,phyASV.f.c.water.M)


############################### Subset by coral species/disturbance ##########################
# Subset by coral species and disturbance level
phyASV.f.c.coral.Peyd <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phyASV.f.c.coral.Peyd <- prune_taxa(taxa_sums(phyASV.f.c.coral.Peyd)>0,phyASV.f.c.coral.Peyd)
phyASV.f.c.coral.MAeq <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Montipora_foliosa")
phyASV.f.c.coral.MAeq <- prune_taxa(taxa_sums(phyASV.f.c.coral.MAeq)>0,phyASV.f.c.coral.MAeq)
phyASV.f.c.coral.Plob <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Porites_lobata")
phyASV.f.c.coral.Plob <- prune_taxa(taxa_sums(phyASV.f.c.coral.Plob)>0,phyASV.f.c.coral.Plob)

phyASV.f.c.coral.Peyd.VH <- subset_samples(phyASV.f.c.coral.Peyd,sample_data(phyASV.f.c.coral.Peyd)$Dist=="VeryHigh")
phyASV.f.c.coral.Peyd.VH <- prune_taxa(taxa_sums(phyASV.f.c.coral.Peyd.VH)>0,phyASV.f.c.coral.Peyd.VH)
phyASV.f.c.coral.MAeq.VH <- subset_samples(phyASV.f.c.coral.MAeq,sample_data(phyASV.f.c.coral.MAeq)$Dist=="VeryHigh")
phyASV.f.c.coral.MAeq.VH <- prune_taxa(taxa_sums(phyASV.f.c.coral.MAeq.VH)>0,phyASV.f.c.coral.MAeq.VH)
phyASV.f.c.coral.Plob.VH <- subset_samples(phyASV.f.c.coral.Plob,sample_data(phyASV.f.c.coral.Plob)$Dist=="VeryHigh")
phyASV.f.c.coral.Plob.VH <- prune_taxa(taxa_sums(phyASV.f.c.coral.Plob.VH)>0,phyASV.f.c.coral.Plob.VH)

phyASV.f.c.coral.Peyd.M <- subset_samples(phyASV.f.c.coral.Peyd,sample_data(phyASV.f.c.coral.Peyd)$Dist=="Medium")
phyASV.f.c.coral.Peyd.M <- prune_taxa(taxa_sums(phyASV.f.c.coral.Peyd.M)>0,phyASV.f.c.coral.Peyd.M)
phyASV.f.c.coral.MAeq.M <- subset_samples(phyASV.f.c.coral.MAeq,sample_data(phyASV.f.c.coral.MAeq)$Dist=="Medium")
phyASV.f.c.coral.MAeq.M <- prune_taxa(taxa_sums(phyASV.f.c.coral.MAeq.M)>0,phyASV.f.c.coral.MAeq.M)
phyASV.f.c.coral.Plob.M <- subset_samples(phyASV.f.c.coral.Plob,sample_data(phyASV.f.c.coral.Plob)$Dist=="Medium")
phyASV.f.c.coral.Plob.M <- prune_taxa(taxa_sums(phyASV.f.c.coral.Plob.M)>0,phyASV.f.c.coral.Plob.M)

############################# Subset by Field Season #################################
# Subset by field season
sediment.before <- subset_samples(phyASV.f.c.sediment, data.frame(sample_data(phyASV.f.c.sediment))$field_season == "KI2014",prune=TRUE)
sediment.before <- subset_taxa(sediment.before, taxa_sums(sediment.before) > 0, prune=TRUE)

sediment.storm <- subset_samples(phyASV.f.c.sediment, data.frame(sample_data(phyASV.f.c.sediment))$field_season == "KI2015a_Post", prune=TRUE)
sediment.storm <- subset_taxa(sediment.storm, taxa_sums(sediment.storm) > 0, prune=TRUE)

sediment.after <- subset_samples(phyASV.f.c.sediment, data.frame(sample_data(phyASV.f.c.sediment))$field_season == "KI2015b", prune=TRUE)
sediment.after <- subset_taxa(sediment.after, taxa_sums(sediment.after) > 0, prune=TRUE)

water.before <- subset_samples(phyASV.f.c.water, data.frame(sample_data(phyASV.f.c.water))$field_season == "KI2014",prune=TRUE)
water.before <- subset_taxa(water.before, taxa_sums(water.before) > 0, prune=TRUE)

water.storm <- subset_samples(phyASV.f.c.water, data.frame(sample_data(phyASV.f.c.water))$field_season == "KI2015a_Post", prune=TRUE)
water.storm <- subset_taxa(water.storm, taxa_sums(water.storm) > 0, prune=TRUE)

water.after <- subset_samples(phyASV.f.c.water, data.frame(sample_data(phyASV.f.c.water))$field_season == "KI2015b", prune=TRUE)
water.after <- subset_taxa(water.after, taxa_sums(water.after) > 0, prune=TRUE)

coral.before <- subset_samples(phyASV.f.c.coral, data.frame(sample_data(phyASV.f.c.coral))$field_season == "KI2014",prune=TRUE)
coral.before <- subset_taxa(coral.before, taxa_sums(coral.before) > 0, prune=TRUE)

coral.storm <- subset_samples(phyASV.f.c.coral, data.frame(sample_data(phyASV.f.c.coral))$field_season == "KI2015a_Post", prune=TRUE)
coral.storm <- subset_taxa(coral.storm, taxa_sums(coral.storm) > 0, prune=TRUE)

coral.after <- subset_samples(phyASV.f.c.coral, data.frame(sample_data(phyASV.f.c.coral))$field_season == "KI2015b", prune=TRUE)
coral.after <- subset_taxa(coral.after, taxa_sums(coral.after) > 0, prune=TRUE)

############################ Types and subclades by compartment #########################
all.types <- unique(data.frame(tax_table(phyASV.f.c))$hit)
sediment.types <- unique(data.frame(tax_table(phyASV.f.c.sediment))$hit)
water.types <- unique(data.frame(tax_table(phyASV.f.c.water))$hit)
coral.types <- unique(data.frame(tax_table(phyASV.f.c.coral))$hit)

sediment.types.subclade <- sediment.types
sediment.types.subclade <- gsub("_.*","",sediment.types.subclade)
sediment.types.subclade <- gsub("\\..*","",sediment.types.subclade)
sediment.types.subclade <- unique(sediment.types.subclade)

water.types.subclade <- water.types
water.types.subclade <- gsub("_.*","",water.types.subclade)
water.types.subclade <- gsub("\\..*","",water.types.subclade)
water.types.subclade <- unique(water.types.subclade)

coral.types.subclade <- coral.types
coral.types.subclade <- gsub("_.*","",coral.types.subclade)
coral.types.subclade <- gsub("\\..*","",coral.types.subclade)
coral.types.subclade <- unique(coral.types.subclade)

# Types and subclades by compartment and time point
sediment.before.types <- unique(data.frame(tax_table(sediment.before))$hit)
sediment.storm.types <- unique(data.frame(tax_table(sediment.storm))$hit)
sediment.after.types <- unique(data.frame(tax_table(sediment.after))$hit)

sediment.before.types.subclade <- sediment.before.types
sediment.before.types.subclade <- gsub("_.*","",sediment.before.types.subclade)
sediment.before.types.subclade <- gsub("\\..*","",sediment.before.types.subclade)
sediment.before.types.subclade <- unique(sediment.before.types.subclade)

sediment.storm.types.subclade <- sediment.storm.types
sediment.storm.types.subclade <- gsub("_.*","",sediment.storm.types.subclade)
sediment.storm.types.subclade <- gsub("\\..*","",sediment.storm.types.subclade)
sediment.storm.types.subclade <- unique(sediment.storm.types.subclade)

sediment.after.types.subclade <- sediment.after.types
sediment.after.types.subclade <- gsub("_.*","",sediment.after.types.subclade)
sediment.after.types.subclade <- gsub("\\..*","",sediment.after.types.subclade)
sediment.after.types.subclade <- unique(sediment.after.types.subclade)


water.before.types <- unique(data.frame(tax_table(water.before))$hit)
water.storm.types <- unique(data.frame(tax_table(water.storm))$hit)
water.after.types <- unique(data.frame(tax_table(water.after))$hit)

water.before.types.subclade <- water.before.types
water.before.types.subclade <- gsub("_.*","",water.before.types.subclade)
water.before.types.subclade <- gsub("\\..*","",water.before.types.subclade)
water.before.types.subclade <- unique(water.before.types.subclade)

water.storm.types.subclade <- water.storm.types
water.storm.types.subclade <- gsub("_.*","",water.storm.types.subclade)
water.storm.types.subclade <- gsub("\\..*","",water.storm.types.subclade)
water.storm.types.subclade <- unique(water.storm.types.subclade)

water.after.types.subclade <- water.after.types
water.after.types.subclade <- gsub("_.*","",water.after.types.subclade)
water.after.types.subclade <- gsub("\\..*","",water.after.types.subclade)
water.after.types.subclade <- unique(water.after.types.subclade)

coral.before.types <- unique(data.frame(tax_table(coral.before))$hit)
coral.storm.types <- unique(data.frame(tax_table(coral.storm))$hit)
coral.after.types <- unique(data.frame(tax_table(coral.after))$hit)

coral.before.types.subclade <- coral.before.types
coral.before.types.subclade <- gsub("_.*","",coral.before.types.subclade)
coral.before.types.subclade <- gsub("\\..*","",coral.before.types.subclade)
coral.before.types.subclade <- unique(coral.before.types.subclade)

coral.storm.types.subclade <- coral.storm.types
coral.storm.types.subclade <- gsub("_.*","",coral.storm.types.subclade)
coral.storm.types.subclade <- gsub("\\..*","",coral.storm.types.subclade)
coral.storm.types.subclade <- unique(coral.storm.types.subclade)

coral.after.types.subclade <- coral.after.types
coral.after.types.subclade <- gsub("_.*","",coral.after.types.subclade)
coral.after.types.subclade <- gsub("\\..*","",coral.after.types.subclade)
coral.after.types.subclade <- unique(coral.after.types.subclade)


####################### Types and subclades by coral species #####################
Peyd.types <- unique(data.frame(tax_table(phyASV.f.c.coral.Peyd))$hit)
Peyd.types.subclade <- Peyd.types
Peyd.types.subclade <- gsub("_.*","",Peyd.types.subclade)
Peyd.types.subclade <- gsub("\\..*","",Peyd.types.subclade)
Peyd.types.subclade <- unique(Peyd.types.subclade)

Plob.types <- unique(data.frame(tax_table(phyASV.f.c.coral.Plob))$hit)
Plob.types.subclade <- Plob.types
Plob.types.subclade <- gsub("_.*","",Plob.types.subclade)
Plob.types.subclade <- gsub("\\..*","",Plob.types.subclade)
Plob.types.subclade <- unique(Plob.types.subclade)

MAeq.types <- unique(data.frame(tax_table(phyASV.f.c.coral.MAeq))$hit)
MAeq.types.subclade <- MAeq.types
MAeq.types.subclade <- gsub("_.*","",MAeq.types.subclade)
MAeq.types.subclade <- gsub("\\..*","",MAeq.types.subclade)
MAeq.types.subclade <- unique(MAeq.types.subclade)

############## Types and subclades by compartment and disturbance level ################
sediment.VH.types <- unique(data.frame(tax_table(phyASV.f.c.sediment.VH))$hit)
sediment.M.types <- unique(data.frame(tax_table(phyASV.f.c.sediment.M))$hit)
water.VH.types <- unique(data.frame(tax_table(phyASV.f.c.water.VH))$hit)
water.M.types <- unique(data.frame(tax_table(phyASV.f.c.water.M))$hit)
coral.VH.types <- unique(data.frame(tax_table(phyASV.f.c.coral.VH))$hit)
coral.M.types <- unique(data.frame(tax_table(phyASV.f.c.coral.M))$hit)

sediment.VH.types.subclade <- sediment.VH.types
sediment.VH.types.subclade <- gsub("_.*","",sediment.VH.types.subclade)
sediment.VH.types.subclade <- gsub("\\..*","",sediment.VH.types.subclade)
sediment.VH.types.subclade <- unique(sediment.VH.types.subclade)
sediment.M.types.subclade <- sediment.M.types
sediment.M.types.subclade <- gsub("_.*","",sediment.M.types.subclade)
sediment.M.types.subclade <- gsub("\\..*","",sediment.M.types.subclade)
sediment.M.types.subclade <- unique(sediment.M.types.subclade)

water.VH.types.subclade <- water.VH.types
water.VH.types.subclade <- gsub("_.*","",water.VH.types.subclade)
water.VH.types.subclade <- gsub("\\..*","",water.VH.types.subclade)
water.VH.types.subclade <- unique(water.VH.types.subclade)
water.M.types.subclade <- water.M.types
water.M.types.subclade <- gsub("_.*","",water.M.types.subclade)
water.M.types.subclade <- gsub("\\..*","",water.M.types.subclade)
water.M.types.subclade <- unique(water.M.types.subclade)

coral.VH.types.subclade <- coral.VH.types
coral.VH.types.subclade <- gsub("_.*","",coral.VH.types.subclade)
coral.VH.types.subclade <- gsub("\\..*","",coral.VH.types.subclade)
coral.VH.types.subclade <- unique(coral.VH.types.subclade)
coral.M.types.subclade <- coral.M.types
coral.M.types.subclade <- gsub("_.*","",coral.M.types.subclade)
coral.M.types.subclade <- gsub("\\..*","",coral.M.types.subclade)
coral.M.types.subclade <- unique(coral.M.types.subclade)

################## Denovo by compartment and coral species ###############
# By denovo
all.denovo <- unique(data.frame(tax_table(phyASV.f.c))$otu)
sediment.denovo <- unique(data.frame(tax_table(phyASV.f.c.sediment))$otu)
water.denovo <- unique(data.frame(tax_table(phyASV.f.c.water))$otu)
coral.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral))$otu)

sediment.denovo.subclade <- sediment.denovo
sediment.denovo.subclade <- unique(sediment.denovo.subclade)

water.denovo.subclade <- water.denovo
water.denovo.subclade <- unique(water.denovo.subclade)

coral.denovo.subclade <- coral.denovo
coral.denovo.subclade <- unique(coral.denovo.subclade)

Peyd.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.Peyd))$otu)
Peyd.denovo.subclade <- Peyd.denovo
Peyd.denovo.subclade <- unique(Peyd.denovo.subclade)

Plob.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.Plob))$otu)
Plob.denovo.subclade <- Plob.denovo
Plob.denovo.subclade <- unique(Plob.denovo.subclade)

MAeq.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.MAeq))$otu)
MAeq.denovo.subclade <- MAeq.denovo
MAeq.denovo.subclade <- unique(MAeq.denovo.subclade)

############# Denovo by compartment/coral species and disturbance ###############
sediment.VH.denovo <- unique(data.frame(tax_table(phyASV.f.c.sediment.VH))$otu)
sediment.M.denovo <- unique(data.frame(tax_table(phyASV.f.c.sediment.M))$otu)
water.VH.denovo <- unique(data.frame(tax_table(phyASV.f.c.water.VH))$otu)
water.M.denovo <- unique(data.frame(tax_table(phyASV.f.c.water.M))$otu)
coral.VH.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.VH))$otu)
coral.M.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.M))$otu)

sediment.VH.denovo.subclade <- unique(sediment.VH.denovo)
sediment.M.denovo.subclade <- unique(sediment.M.denovo)
water.VH.denovo.subclade <- unique(water.VH.denovo)
water.M.denovo.subclade <- unique(water.M.denovo)
coral.VH.denovo.subclade <- unique(coral.VH.denovo)
coral.M.denovo.subclade <- unique(coral.M.denovo)

# denovo and subclades by compartment and disturbance level
Peyd.VH.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.Peyd.VH))$hit)
Peyd.M.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.Peyd.M))$hit)
Plob.VH.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.Plob.VH))$hit)
Plob.M.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.Plob.M))$hit)
MAeq.VH.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.MAeq.VH))$hit)
MAeq.M.denovo <- unique(data.frame(tax_table(phyASV.f.c.coral.MAeq.M))$hit)

Peyd.VH.denovo.subclade <- unique(Peyd.VH.denovo)
Peyd.M.denovo.subclade <- unique(Peyd.M.denovo)
Plob.VH.denovo.subclade <- unique(Plob.VH.denovo)
Plob.M.denovo.subclade <- unique(Plob.M.denovo)
MAeq.VH.denovo.subclade <- unique(MAeq.VH.denovo)
MAeq.M.denovo.subclade <- unique(MAeq.M.denovo)

#################### Save grouped data as RData file ##########################
save(list=ls(),file="data/KI_Compartment_f_coral_grouped.RData")
save(list=c("phyASV.f.c.p"),file="data/KI_Compartment_phyASVfcp.RData")
