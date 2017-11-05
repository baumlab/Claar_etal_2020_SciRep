# Import Libraries
library(stringr)
library(reshape2)
library(phyloseq)
library(seqinr)
library(phangorn) 
library(caroline)

# Clear workspace
rm(list=ls())

# Load filtered RData object from output of filter_notsym.R script
load("analyses/KI_Compartment.RData")

# Filter OTUs by minimum count
# Set threshold count
n <- 5
# Identify OTUs below threshold count
taxa <- taxa_sums(phy.f)[which(taxa_sums(phy.f) >= n)]
# Remove taxa below threshold count
phy.f <- prune_taxa(names(taxa), phy.f)

# Filter samples by minimum count
# Set threshold number of reads
sn <- 200
# Remove samples with fewer reads than threshold
phy.f <- prune_samples(sample_sums(phy.f)>=sn, phy.f)

# Characterize sites by disturbance level
VeryHigh <- c(33,40,32,31,27,30,26)
High <- c(25,3,38,24)
HighMed <- c(9,34,35,8,14,6)
LowMed <- c(2,22,1,23)
Low <- c(7,13,12,4,36,5,37)
VeryLow <- c(10,21,11,20,16,15,39,19,18,17)

sample_data(phy.f)$site <- factor(sample_data(phy.f)$site)
sample_data(phy.f)$Dist <- sample_data(phy.f)$site

for (i in VeryHigh){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryHigh"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("VeryHigh")))
}

for (i in High){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("High"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("High")))
}

for (i in HighMed){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("HighMed"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("HighMed")))
}

for (i in LowMed){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("LowMed"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("LowMed")))
}

for (i in Low){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Low"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("Low")))
}

for (i in VeryLow){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryLow"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("VeryLow")))
}

levels(sample_data(phy.f)$Dist) <- c("VeryHigh","High","HighMed","Low","VeryLow")

# Assign new name for clarity
phy97.f.c <- phy.f

# Make a tax_table column for "clade"
# Rename rank_names (aka column names in phy97.f.c)
colnames(tax_table(phy97.f.c)) <- c(otu = "otu", sim = "sim", del = "del", ins = "ins", mis = "mis", len = "len", score = "score", hit = "hit", i = "clade")
# Copy "hit" into last column
tax_table(phy97.f.c)[,9] <- tax_table(phy97.f.c)[,8]

# Use gsub to replace full "hit" name with only clade letter
tax_table(phy97.f.c)[,9] <- gsub("^A.*", "A", tax_table(phy97.f.c)[,9])
tax_table(phy97.f.c)[,9] <- gsub("^B.*", "B", tax_table(phy97.f.c)[,9])
tax_table(phy97.f.c)[,9] <- gsub("^C.*", "C", tax_table(phy97.f.c)[,9])
tax_table(phy97.f.c)[,9] <- gsub("^D.*", "D", tax_table(phy97.f.c)[,9])
tax_table(phy97.f.c)[,9] <- gsub("^F.*", "F", tax_table(phy97.f.c)[,9])
tax_table(phy97.f.c)[,9] <- gsub("^G.*", "G", tax_table(phy97.f.c)[,9])
tax_table(phy97.f.c)[,9] <- gsub("^I.*", "I", tax_table(phy97.f.c)[,9])

# Write tax table for phylogenetically-informed diversity
write.table(data.frame(tax_table(phy97.f.c)), "data/tax_table.txt", row.names=T, quote=F)
# Write otu table for phylogenetically-informed diversity
write.delim(data.frame(otu_table(phy97.f.c)), "data/otu_table.tsv", quote = FALSE, row.names = T, sep = "\t")

#https://rdrr.io/rforge/seqinr/man/dist.alignment.html
#returns sqrt of pairwise genetic distance, then squared the matrices
A.seqs <- read.alignment(file = "data/Bioinf/tree/A_tree_seqs_aligned_clean.fasta", format= "fasta")
A.dis <- (as.matrix(dist.alignment(A.seqs, matrix = "identity" )))^2
write.csv(A.dis, file="data/Bioinf/tree/A.dis.matx.csv")

B.seqs <- read.alignment(file = "data/Bioinf/tree/B_tree_seqs_aligned_clean.fasta", format= "fasta")
B.dis <- (as.matrix(dist.alignment(B.seqs, matrix = "identity" )))^2
write.csv(B.dis, file="data/Bioinf/tree/B.dis.matx.csv")

C.seqs <- read.alignment(file = "data/Bioinf/tree/C_tree_seqs_aligned_clean.fasta", format= "fasta")
C.dis <- (as.matrix(dist.alignment(C.seqs, matrix = "identity" )))^2
write.csv(C.dis, file="data/Bioinf/tree/C.dis.matx.csv")

D.seqs <- read.alignment(file = "data/Bioinf/tree/D_tree_seqs_aligned_clean.fasta", format= "fasta")
D.dis <- (as.matrix(dist.alignment(D.seqs, matrix = "identity" )))^2
write.csv(D.dis, file="data/Bioinf/tree/D.dis.matx.csv")

F.seqs <- read.alignment(file = "data/Bioinf/tree/F_tree_seqs_aligned_clean.fasta", format= "fasta")
F.dis <- (as.matrix(dist.alignment(F.seqs, matrix = "identity" )))^2
write.csv(F.dis, file="data/Bioinf/tree/F.dis.matx.csv")

G.seqs <- read.alignment(file = "data/Bioinf/tree/G_tree_seqs_aligned_clean.fasta", format= "fasta")
G.dis <- (as.matrix(dist.alignment(G.seqs, matrix = "identity" )))^2
write.csv(G.dis, file="data/Bioinf/tree/G.dis.matx.csv")

I.seqs <- read.alignment(file = "data/Bioinf/tree/I_tree_seqs_aligned_clean.fasta", format= "fasta")
I.dis <- (as.matrix(dist.alignment(I.seqs, matrix = "identity" )))^2
write.csv(I.dis, file="data/Bioinf/tree/I.dis.matx.csv")

#give clade distances using average 28s distance from Pochon and Gates 2010
A_B <- matrix(0.219, ncol=ncol(A.dis), nrow=nrow(B.dis), dimnames=list(rownames(B.dis), colnames(A.dis)))
A_C <- matrix(0.1960, ncol=ncol(A.dis), nrow=nrow(C.dis), dimnames=list(rownames(C.dis), colnames(A.dis)))
A_D <- matrix(0.1775, ncol=ncol(A.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(A.dis)))
A_F <- matrix(0.2085, ncol=ncol(A.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(A.dis)))
A_G <- matrix(0.216, ncol=ncol(A.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(A.dis)))
A_I <- matrix(0.205, ncol=ncol(A.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(A.dis)))
B_C <- matrix(0.114, ncol=ncol(B.dis), nrow=nrow(C.dis), dimnames=list(rownames(C.dis), colnames(B.dis)))
B_D <- matrix(0.1705, ncol=ncol(B.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(B.dis)))
B_F <- matrix(0.1355, ncol=ncol(B.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(B.dis)))
B_G <- matrix(0.21, ncol=ncol(B.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(B.dis)))
B_I <- matrix(0.151, ncol=ncol(B.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(B.dis)))
C_D <- matrix(0.1520, ncol=ncol(C.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(C.dis)))
C_F <- matrix(0.0945, ncol=ncol(C.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(C.dis)))
C_G <- matrix(0.187, ncol=ncol(C.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(C.dis)))
C_I <- matrix(0.137, ncol=ncol(C.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(C.dis)))
D_F <- matrix(0.1691, ncol=ncol(D.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(D.dis)))
D_G <- matrix(0.1795, ncol=ncol(D.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(D.dis)))
D_I <- matrix(0.169125, ncol=ncol(D.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(D.dis)))
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

phy_tree(uber.tree)

# Slot uber tree into the phy_tree slot of the phyloseq object
phy_tree(phy97.f.c) <- phy_tree(uber.tree)

# Transform sample counts to proportional abundance for downstream analyses
phy97.f.c.p <- transform_sample_counts(phy97.f.c, function(x) x/sum(x))

# Note, this is now the same as phy97.f.c - renaming and subsetting was completed before when water and sediment samples were included
phy97.f.c.coral <- subset_samples(phy97.f.c,SampleType=="coral")
phy97.f.c.coral <- subset_taxa(phy97.f.c.coral, taxa_sums(phy97.f.c.coral) > 0, prune=TRUE)

phy97.f.c.water <- subset_samples(phy97.f.c,SampleType=="water")
phy97.f.c.water <- subset_taxa(phy97.f.c.water, taxa_sums(phy97.f.c.water) > 0, prune=TRUE)

phy97.f.c.sediment <- subset_samples(phy97.f.c,SampleType=="sediment")
phy97.f.c.sediment <- subset_taxa(phy97.f.c.sediment, taxa_sums(phy97.f.c.sediment) > 0, prune=TRUE)

# Transform sample counts to proportional abundance for downstream analyses
phy97.f.c.coral.p <- transform_sample_counts(phy97.f.c.coral, function(x) x/sum(x))
phy97.f.c.water.p <- transform_sample_counts(phy97.f.c.water, function(x) x/sum(x))
phy97.f.c.sediment.p <- transform_sample_counts(phy97.f.c.sediment, function(x) x/sum(x))

# Calculate number of sequences in phy97.f.c
total_seqs <- sum(taxa_sums(phy97.f.c))
coral_seqs <- sum(taxa_sums(phy97.f.c.coral))
water_seqs <- sum(taxa_sums(phy97.f.c.water))
sediment_seqs <- sum(taxa_sums(phy97.f.c.sediment))

# Save grouped data as RData file
save(list=ls(pattern="phy97.f.c."), sediment_seqs, total_seqs, coral_seqs, water_seqs, file = "data/KI_Compartment_f_coral_grouped.RData")
