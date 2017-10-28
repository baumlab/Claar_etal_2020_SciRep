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
load("analyses/KI_Platy.RData")

phy.f.coral <- phy.f

# Filter OTUs by minimum count
# Set threshold count
n <- 5
# Identify OTUs below threshold count
taxa <- taxa_sums(phy.f.coral)[which(taxa_sums(phy.f.coral) >= n)]
# Remove taxa below threshold count
phy.f.coral <- prune_taxa(names(taxa), phy.f.coral)

# Filter samples by minimum count
# Set threshold number of reads
sn <- 200
# Remove samples with fewer reads than threshold
phy.f.coral <- prune_samples(sample_sums(phy.f.coral)>=sn, phy.f.coral)

# I decided not to do this, because if they were present in > the threshold before, they are most likely 'real' otus. Although the samples that they were in didn't sequence very well, that doesn't necessarily mean that they are incorrect. 
# # Filter OTUs by minimum count again in case any dropped below threshold after filtering samples
# # Identify OTUs below threshold count
# taxa <- taxa_sums(phy.f.coral)[which(taxa_sums(phy.f.coral) >= n)]
# # Remove taxa below threshold count
# phy.f.coral <- prune_taxa(names(taxa), phy.f.coral)

# Characterize sites by disturbance level
VeryHigh <- c(33,40,32,31,27,30,26)
High <- c(25,3,38,24)
HighMed <- c(9,34,35,8,14,6)
LowMed <- c(2,22,1,23)
Low <- c(7,13,12,4,36,5,37)
VeryLow <- c(10,21,11,20,16,15,39,19,18,17)

sample_data(phy.f.coral)$site <- factor(sample_data(phy.f.coral)$site)
sample_data(phy.f.coral)$Dist <- sample_data(phy.f.coral)$site

for (i in VeryHigh){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f.coral))$site)))
    (sample_data(phy.f.coral)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryHigh"), as.character(data.frame(sample_data(phy.f.coral))$Dist), as.character("VeryHigh")))
}

for (i in High){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f.coral))$site)))
    (sample_data(phy.f.coral)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("High"), as.character(data.frame(sample_data(phy.f.coral))$Dist), as.character("High")))
}

for (i in HighMed){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f.coral))$site)))
    (sample_data(phy.f.coral)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("HighMed"), as.character(data.frame(sample_data(phy.f.coral))$Dist), as.character("HighMed")))
}

for (i in LowMed){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f.coral))$site)))
    (sample_data(phy.f.coral)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("LowMed"), as.character(data.frame(sample_data(phy.f.coral))$Dist), as.character("LowMed")))
}

for (i in Low){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f.coral))$site)))
    (sample_data(phy.f.coral)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Low"), as.character(data.frame(sample_data(phy.f.coral))$Dist), as.character("Low")))
}

for (i in VeryLow){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f.coral))$site)))
    (sample_data(phy.f.coral)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryLow"), as.character(data.frame(sample_data(phy.f.coral))$Dist), as.character("VeryLow")))
}

levels(sample_data(phy.f.coral)$Dist) <- c("VeryHigh","High","HighMed","Low","VeryLow")

# Assign new name for clarity
phy97.f.c <- phy.f.coral

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

C.seqs <- read.alignment(file = "data/Bioinf/tree/C_tree_seqs_aligned_clean.fasta", format= "fasta")
C.dis <- (as.matrix(dist.alignment(C.seqs, matrix = "identity" )))^2
write.csv(C.dis, file="data/Bioinf/tree/C.dis.matx.csv")

D.seqs <- read.alignment(file = "data/Bioinf/tree/D_tree_seqs_aligned_clean.fasta", format= "fasta")
D.dis <- (as.matrix(dist.alignment(D.seqs, matrix = "identity" )))^2
write.csv(D.dis, file="data/Bioinf/tree/D.dis.matx.csv")

# F.seqs <- read.alignment(file = "data/Bioinf/tree/F_tree_seqs_aligned_clean.fasta", format= "fasta")
# F.dis <- (as.matrix(dist.alignment(F.seqs, matrix = "identity" )))^2
# write.csv(F.dis, file="data/Bioinf/tree/F.dis.matx.csv")

G.seqs <- read.alignment(file = "data/Bioinf/tree/G_tree_seqs_aligned_clean.fasta", format= "fasta")
G.dis <- (as.matrix(dist.alignment(G.seqs, matrix = "identity" )))^2
write.csv(G.dis, file="data/Bioinf/tree/G.dis.matx.csv")

#give clade distances using average 28s distance from Pochon and Gates 2010
A_C <- matrix(0.1960, ncol=ncol(A.dis), nrow=nrow(C.dis), dimnames=list(rownames(C.dis), colnames(A.dis)))
A_D <- matrix(0.1775, ncol=ncol(A.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(A.dis)))
# A_F <- matrix(0.2085, ncol=ncol(A.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(A.dis)))
A_G <- matrix(0.216, ncol=ncol(A.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(A.dis)))
C_D <- matrix(0.1520, ncol=ncol(C.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(C.dis)))
# C_F <- matrix(0.0945, ncol=ncol(C.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(C.dis)))
C_G <- matrix(0.187, ncol=ncol(C.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(C.dis)))
# D_F <- matrix(0.1691, ncol=ncol(D.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(D.dis)))
D_G <- matrix(0.1795, ncol=ncol(D.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(D.dis)))
# F_G <- matrix(0.2072, ncol=ncol(F.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(F.dis)))


#build ACDG matrix
col1 <- rbind(A.dis, A_C, A_D, A_G)
col2 <- rbind(matrix(NA, nrow=nrow(A.dis), ncol=ncol(C.dis), dimnames=list(rownames(A.dis), colnames(C.dis))), C.dis, C_D, C_G)
col3 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(C.dis), ncol=ncol(D.dis), dimnames=list(c(rownames(A.dis), rownames(C.dis)), colnames(D.dis))), D.dis, D_G)
col4 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(C.dis)+nrow(D.dis), ncol=ncol(G.dis), dimnames=list(c(rownames(A.dis), rownames(C.dis), rownames(D.dis)), colnames(G.dis))), G.dis)

ubermatrix <- cbind(col1, col2, col3, col4)
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

# Subset coral data set into individual genus data sets
phy97.f.c.platy <- subset_samples(phy97.f.c,Coral_Species=="Platygyra_sp")
phy97.f.c.platy <- subset_taxa(phy97.f.c.platy, taxa_sums(phy97.f.c.platy) > 0, prune=TRUE)
phy97.f.c.platy.AD <- subset_samples(phy97.f.c.platy,Status=="alive"|Status=="dead")

phy97.f.c.fpenta <- subset_samples(phy97.f.c,Coral_Species=="Favites_pentagona")
phy97.f.c.fpenta <- subset_taxa(phy97.f.c.fpenta, taxa_sums(phy97.f.c.fpenta) > 0, prune=TRUE)
phy97.f.c.fpenta.AD <- subset_samples(phy97.f.c.fpenta,Status=="alive"|Status=="dead")

# Note, this is now the same as phy97.f.c - renaming and subsetting was completed before when water and sediment samples were included
phy97.f.c.coral <- subset_samples(phy97.f.c,SampleType=="coral")
phy97.f.c.coral <- subset_taxa(phy97.f.c.coral, taxa_sums(phy97.f.c.coral) > 0, prune=TRUE)

# Subset coral samples that have been confirmed either alive or dead (remove unknowns/gone)
phy97.f.c.coral.AD <- subset_samples(phy97.f.c.coral,Status=="alive"|Status=="dead")

# Transform sample counts to proportional abundance for downstream analyses
phy97.f.c.platy.p <- transform_sample_counts(phy97.f.c.platy, function(x) x/sum(x))
# Transform sample counts to proportional abundance for downstream analyses
phy97.f.c.fpenta.p <- transform_sample_counts(phy97.f.c.fpenta, function(x) x/sum(x))

# Subset coral by site
for (i in unique(data.frame(sample_data(phy97.f.c.coral))$site)){
  nam <- paste("phy97.f.c.coral.site",i,sep="")
  a <- eval(subset_samples(phy97.f.c.coral,site==i, prune=TRUE))
  assign(nam,a)
  b <- eval(subset_taxa(a, taxa_sums(a) > 0, prune=TRUE))
  assign(nam,b)
}

# Subset platygyra only by site
for (i in unique(data.frame(sample_data(phy97.f.c.platy))$site)){
  nam <- paste("phy97.f.c.platy.site",i,sep="")
  a <- eval(subset_samples(phy97.f.c.platy,site==i, prune=TRUE))
  assign(nam,a)
  b <- eval(subset_taxa(a, taxa_sums(a) > 0, prune=TRUE))
  assign(nam,b)
}

# Subset fpenta only by site
for (i in unique(data.frame(sample_data(phy97.f.c.fpenta))$site)){
  nam <- paste("phy97.f.c.fpenta.site",i,sep="")
  a <- eval(subset_samples(phy97.f.c.fpenta,site==i, prune=TRUE))
  assign(nam,a)
  b <- eval(subset_taxa(a, taxa_sums(a) > 0, prune=TRUE))
  assign(nam,b)
}

# Subset coral by disturbance level
for (i in unique(data.frame(sample_data(phy97.f.c.coral))$Dist)){
  nam <- paste("phy97.f.c.coral.",i,sep="")
  a <- eval(subset_samples(phy97.f.c.coral,Dist==i, prune=TRUE))
  assign(nam,a)
  b <- eval(subset_taxa(a, taxa_sums(a) > 0, prune=TRUE))
  assign(nam,b)
}

# Subset platygyra only by Disturbance level
for (i in unique(data.frame(sample_data(phy97.f.c.platy))$Dist)){
  nam <- paste("phy97.f.c.platy.",i,sep="")
  a <- eval(subset_samples(phy97.f.c.platy,Dist==i, prune=TRUE))
  assign(nam,a)
  b <- eval(subset_taxa(a, taxa_sums(a) > 0, prune=TRUE))
  assign(nam,b)
}

# Subset fpenta only by disturbance level
for (i in unique(data.frame(sample_data(phy97.f.c.fpenta))$Dist)){
  nam <- paste("phy97.f.c.fpenta.",i,sep="")
  a <- eval(subset_samples(phy97.f.c.fpenta,Dist==i, prune=TRUE))
  assign(nam,a)
  b <- eval(subset_taxa(a, taxa_sums(a) > 0, prune=TRUE))
  assign(nam,b)
}


# Transform sample counts to proportional abundance for downstream analyses
for (i in unique(data.frame(sample_data(phy97.f.c.coral))$site)){
  nam <- paste("phy97.f.c.coral.site",i,".p",sep="")
  c <- paste("phy97.f.c.coral.site",i,sep="")
  a <- eval(transform_sample_counts(get(c), function(x) x/sum(x)))
  assign(nam,a)
}

# Subset by Field Season
phy97.f.c.coral.KI2014 <- subset_samples(phy97.f.c.coral,field_season=="KI2014")
phy97.f.c.coral.KI2014 <- subset_taxa(phy97.f.c.coral.KI2014, taxa_sums(phy97.f.c.coral.KI2014) > 0, prune=TRUE)
phy97.f.c.coral.KI2014.AD <- subset_samples(phy97.f.c.coral.KI2014,Status=="alive"|Status=="dead")
phy97.f.c.coral.KI2015a <- subset_samples(phy97.f.c.coral,field_season=="KI2015a")
phy97.f.c.coral.KI2015a <- subset_taxa(phy97.f.c.coral.KI2015a, taxa_sums(phy97.f.c.coral.KI2015a) > 0, prune=TRUE)
phy97.f.c.coral.KI2015a.AD <- subset_samples(phy97.f.c.coral.KI2015a,Status=="alive"|Status=="dead")
phy97.f.c.coral.KI2015b <- subset_samples(phy97.f.c.coral,field_season=="KI2015b")
phy97.f.c.coral.KI2015b <- subset_taxa(phy97.f.c.coral.KI2015b, taxa_sums(phy97.f.c.coral.KI2015b) > 0, prune=TRUE)
phy97.f.c.coral.KI2015b.AD <- subset_samples(phy97.f.c.coral.KI2015b,Status=="alive"|Status=="dead")
phy97.f.c.coral.KI2015c <- subset_samples(phy97.f.c.coral,field_season=="KI2015c")
phy97.f.c.coral.KI2015c <- subset_taxa(phy97.f.c.coral.KI2015c, taxa_sums(phy97.f.c.coral.KI2015c) > 0, prune=TRUE)
phy97.f.c.coral.KI2015c.AD <- subset_samples(phy97.f.c.coral.KI2015c,Status=="alive"|Status=="dead")
phy97.f.c.coral.KI2016a <- subset_samples(phy97.f.c.coral,field_season=="KI2016a")
phy97.f.c.coral.KI2016a <- subset_taxa(phy97.f.c.coral.KI2016a, taxa_sums(phy97.f.c.coral.KI2016a) > 0, prune=TRUE)
phy97.f.c.coral.KI2016a.AD <- subset_samples(phy97.f.c.coral.KI2016a,Status=="alive"|Status=="dead")

# Subset coral samples to only keep samples taken before the event (==KI2014 to May 2015)
phy97.f.c.coral.AD.before <- subset_samples(phy97.f.c.coral.AD,field_season!="KI2016a", prune=TRUE)
phy97.f.c.coral.AD.before <- subset_samples(phy97.f.c.coral.AD.before,field_season!="KI2015c", prune=TRUE)
# Subset Platygyra samples to only keep samples taken before the event (==KI2014 to May 2015)
phy97.f.c.platy.AD.before <- subset_samples(phy97.f.c.platy.AD,field_season!="KI2016a", prune=TRUE)
phy97.f.c.platy.AD.before <- subset_samples(phy97.f.c.platy.AD.before,field_season!="KI2015c", prune=TRUE)
# Subset fpenta samples to only keep samples taken before the event (==KI2014 to May 2015)
phy97.f.c.fpenta.AD.before <- subset_samples(phy97.f.c.fpenta.AD,field_season!="KI2016a", prune=TRUE)
phy97.f.c.fpenta.AD.before <- subset_samples(phy97.f.c.fpenta.AD.before,field_season!="KI2015c", prune=TRUE)
# Subset coral samples to only keep samples taken during/after the event (== July 2015 - 2016)
phy97.f.c.coral.AD.da <- subset_samples(phy97.f.c.coral.AD,field_season=="KI2016a"|field_season=="KI2015c", prune=TRUE)
# Subset Platygyra samples to only keep samples taken during/after the event (== July 2015 - 2016)
phy97.f.c.platy.AD.da <- subset_samples(phy97.f.c.platy.AD,field_season=="KI2016a"|field_season=="KI2015c", prune=TRUE)
# Subset FPenta samples to only keep samples taken during/after the event (== July 2015 - 2016)
phy97.f.c.fpenta.AD.da <- subset_samples(phy97.f.c.fpenta.AD,field_season=="KI2016a"|field_season=="KI2015c", prune=TRUE)
# Subset coral samples to only keep samples taken after the event (== 2016)
phy97.f.c.coral.AD.after <- subset_samples(phy97.f.c.coral.AD,field_season=="KI2016a", prune=TRUE)
# Subset Platygyra samples to only keep samples taken after the event (== 2016)
phy97.f.c.platy.AD.after <- subset_samples(phy97.f.c.platy.AD,field_season=="KI2016a", prune=TRUE)
# Subset FPenta samples to only keep samples taken during/after the event (== 2016)
phy97.f.c.fpenta.AD.after <- subset_samples(phy97.f.c.fpenta.AD,field_season=="KI2016a", prune=TRUE)
# Subset coral samples to only keep samples taken during the event (== July 2015)
phy97.f.c.coral.AD.during <- subset_samples(phy97.f.c.coral.AD,field_season=="KI2015c", prune=TRUE)
# Subset Platygyra samples to only keep samples taken during the event (== July 2015)
phy97.f.c.platy.AD.during <- subset_samples(phy97.f.c.platy.AD,field_season=="KI2015c", prune=TRUE)
# Subset FPenta samples to only keep samples taken during the event (== July 2015)
phy97.f.c.fpenta.AD.during <- subset_samples(phy97.f.c.fpenta.AD,field_season=="KI2015c", prune=TRUE)

phy97.f.c.fpenta.AD.before.LVL <- subset_samples(phy97.f.c.fpenta.AD.before, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.fpenta.AD.before.LVL <- subset_taxa(phy97.f.c.fpenta.AD.before.LVL, taxa_sums(phy97.f.c.fpenta.AD.before.LVL) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.before.HM <- subset_samples(phy97.f.c.fpenta.AD.before, Dist=="HighMed", prune=TRUE)
phy97.f.c.fpenta.AD.before.HM <- subset_taxa(phy97.f.c.fpenta.AD.before.HM, taxa_sums(phy97.f.c.fpenta.AD.before.HM) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.before.HVH <- subset_samples(phy97.f.c.fpenta.AD.before, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
phy97.f.c.fpenta.AD.before.LVL <- subset_taxa(phy97.f.c.fpenta.AD.before.LVL, taxa_sums(phy97.f.c.fpenta.AD.before.LVL) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.during.LVL <- subset_samples(phy97.f.c.fpenta.AD.during, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.fpenta.AD.before.LVL <- subset_taxa(phy97.f.c.fpenta.AD.before.LVL, taxa_sums(phy97.f.c.fpenta.AD.before.LVL) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.during.HM <- subset_samples(phy97.f.c.fpenta.AD.during, Dist=="HighMed", prune=TRUE)
phy97.f.c.fpenta.AD.during.HM <- subset_taxa(phy97.f.c.fpenta.AD.during.HM, taxa_sums(phy97.f.c.fpenta.AD.during.HM) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.during.HVH <- subset_samples(phy97.f.c.fpenta.AD.during, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
phy97.f.c.fpenta.AD.during.HVH <- subset_taxa(phy97.f.c.fpenta.AD.during.HVH, taxa_sums(phy97.f.c.fpenta.AD.during.HVH) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.after.LVL <- subset_samples(phy97.f.c.fpenta.AD.after, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.fpenta.AD.after.LVL <- subset_taxa(phy97.f.c.fpenta.AD.after.LVL, taxa_sums(phy97.f.c.fpenta.AD.after.LVL) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.after.HM <- subset_samples(phy97.f.c.fpenta.AD.after, Dist=="HighMed", prune=TRUE)
phy97.f.c.fpenta.AD.after.HM <- subset_taxa(phy97.f.c.fpenta.AD.after.HM, taxa_sums(phy97.f.c.fpenta.AD.after.HM) > 0, prune=TRUE)
phy97.f.c.fpenta.AD.after.HVH <- subset_samples(phy97.f.c.fpenta.AD.after, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
phy97.f.c.fpenta.AD.after.HVH <- subset_taxa(phy97.f.c.fpenta.AD.after.HVH, taxa_sums(phy97.f.c.fpenta.AD.after.HVH) > 0, prune=TRUE)

phy97.f.c.platy.AD.before.LVL <- subset_samples(phy97.f.c.platy.AD.before, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.platy.AD.before.LVL <- subset_taxa(phy97.f.c.platy.AD.before.LVL, taxa_sums(phy97.f.c.platy.AD.before.LVL) > 0, prune=TRUE)
phy97.f.c.platy.AD.before.HM <- subset_samples(phy97.f.c.platy.AD.before, Dist=="HighMed", prune=TRUE)
phy97.f.c.platy.AD.before.HM <- subset_taxa(phy97.f.c.platy.AD.before.HM, taxa_sums(phy97.f.c.platy.AD.before.HM) > 0, prune=TRUE)
phy97.f.c.platy.AD.before.HVH <- subset_samples(phy97.f.c.platy.AD.before, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
phy97.f.c.platy.AD.before.HVH <- subset_taxa(phy97.f.c.platy.AD.before.HVH, taxa_sums(phy97.f.c.platy.AD.before.HVH) > 0, prune=TRUE)
phy97.f.c.platy.AD.during.LVL <- subset_samples(phy97.f.c.platy.AD.during, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.platy.AD.during.LVL <- subset_taxa(phy97.f.c.platy.AD.during.LVL, taxa_sums(phy97.f.c.platy.AD.during.LVL) > 0, prune=TRUE)
phy97.f.c.platy.AD.during.HM <- subset_samples(phy97.f.c.platy.AD.during, Dist=="HighMed", prune=TRUE)
phy97.f.c.platy.AD.during.HM <- subset_taxa(phy97.f.c.platy.AD.during.HM, taxa_sums(phy97.f.c.platy.AD.during.HM) > 0, prune=TRUE)
phy97.f.c.platy.AD.during.HVH <- subset_samples(phy97.f.c.platy.AD.during, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
phy97.f.c.platy.AD.during.HVH <- subset_taxa(phy97.f.c.platy.AD.during.HVH, taxa_sums(phy97.f.c.platy.AD.during.HVH) > 0, prune=TRUE)
phy97.f.c.platy.AD.after.LVL <- subset_samples(phy97.f.c.platy.AD.after, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.platy.AD.after.LVL <- subset_taxa(phy97.f.c.platy.AD.after.LVL, taxa_sums(phy97.f.c.platy.AD.after.LVL) > 0, prune=TRUE)
phy97.f.c.platy.AD.after.HM <- subset_samples(phy97.f.c.platy.AD.after, Dist=="HighMed", prune=TRUE)
phy97.f.c.platy.AD.after.HM <- subset_taxa(phy97.f.c.platy.AD.after.HM, taxa_sums(phy97.f.c.platy.AD.after.HM) > 0, prune=TRUE)
phy97.f.c.platy.AD.after.HVH <- subset_samples(phy97.f.c.platy.AD.after, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
phy97.f.c.platy.AD.after.HVH <- subset_taxa(phy97.f.c.platy.AD.after.HVH, taxa_sums(phy97.f.c.platy.AD.after.HVH) > 0, prune=TRUE)

# phy97.f.c.coral.AD.before.p <- transform_sample_counts(phy97.f.c.coral.AD.before, function(x) x/sum(x))
# phy97.f.c.coral.AD.during.p <- transform_sample_counts(phy97.f.c.coral.AD.during, function(x) x/sum(x))
# phy97.f.c.coral.AD.after.p <- transform_sample_counts(phy97.f.c.coral.AD.after, function(x) x/sum(x))
# phy97.f.c.coral.p <- transform_sample_counts(phy97.f.c.coral, function(x) x/sum(x))
phy97.f.c.coral.before <- subset_samples(phy97.f.c.coral, field_season=="KI2014"|field_season=="KI2015a"|field_season=="KI2015b", prune=TRUE)
phy97.f.c.coral.before.LVL <- subset_samples(phy97.f.c.coral.before, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.coral.before.HM <- subset_samples(phy97.f.c.coral.before, Dist=="HighMed", prune=TRUE)
phy97.f.c.coral.before.HVH <- subset_samples(phy97.f.c.coral.before, Dist=="High" | Dist=="VeryHigh", prune=TRUE)

phy97.f.c.coral.during <- subset_samples(phy97.f.c.coral, field_season=="KI2015c", prune=TRUE)
phy97.f.c.coral.during.LVL <- subset_samples(phy97.f.c.coral.during, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.coral.during.HM <- subset_samples(phy97.f.c.coral.during, Dist=="HighMed", prune=TRUE)
phy97.f.c.coral.during.HVH <- subset_samples(phy97.f.c.coral.during, Dist=="High" | Dist=="VeryHigh", prune=TRUE)

phy97.f.c.coral.after <- subset_samples(phy97.f.c.coral, field_season=="KI2016a", prune=TRUE)
phy97.f.c.coral.after.LVL <- subset_samples(phy97.f.c.coral.after, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
phy97.f.c.coral.after.HM <- subset_samples(phy97.f.c.coral.after, Dist=="HighMed", prune=TRUE)
phy97.f.c.coral.after.HVH <- subset_samples(phy97.f.c.coral.after, Dist=="High" | Dist=="VeryHigh", prune=TRUE)

# phy97.f.c.coral.AD.before.LVL <- subset_taxa(phy97.f.c.coral.AD.before.LVL, taxa_sums(phy97.f.c.coral.AD.before.LVL) > 0, prune=TRUE)
# phy97.f.c.coral.AD.before.HM <- subset_samples(phy97.f.c.coral.AD.before, Dist=="HighMed", prune=TRUE)
# phy97.f.c.coral.AD.before.HM <- subset_taxa(phy97.f.c.coral.AD.before.HM, taxa_sums(phy97.f.c.coral.AD.before.HM) > 0, prune=TRUE)
# phy97.f.c.coral.AD.before.HVH <- subset_samples(phy97.f.c.coral.AD.before, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
# phy97.f.c.coral.AD.before.HVH <- subset_taxa(phy97.f.c.coral.AD.before.HVH, taxa_sums(phy97.f.c.coral.AD.before.HVH) > 0, prune=TRUE)
# phy97.f.c.coral.AD.during.LVL <- subset_samples(phy97.f.c.coral.AD.during, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
# phy97.f.c.coral.AD.during.LVL <- subset_taxa(phy97.f.c.coral.AD.during.LVL, taxa_sums(phy97.f.c.coral.AD.during.LVL) > 0, prune=TRUE)
# phy97.f.c.coral.AD.during.HM <- subset_samples(phy97.f.c.coral.AD.during, Dist=="HighMed", prune=TRUE)
# phy97.f.c.coral.AD.during.HM <- subset_taxa(phy97.f.c.coral.AD.during.HM, taxa_sums(phy97.f.c.coral.AD.during.HM) > 0, prune=TRUE)
# phy97.f.c.coral.AD.during.HVH <- subset_samples(phy97.f.c.coral.AD.during, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
# phy97.f.c.coral.AD.during.HVH <- subset_taxa(phy97.f.c.coral.AD.during.HVH, taxa_sums(phy97.f.c.coral.AD.during.HVH) > 0, prune=TRUE)
# phy97.f.c.coral.AD.after.LVL <- subset_samples(phy97.f.c.coral.AD.after, Dist=="Low" | Dist=="VeryLow", prune=TRUE)
# phy97.f.c.coral.AD.after.LVL <- subset_taxa(phy97.f.c.coral.AD.after.LVL, taxa_sums(phy97.f.c.coral.AD.after.LVL) > 0, prune=TRUE)
# phy97.f.c.coral.AD.after.HM <- subset_samples(phy97.f.c.coral.AD.after, Dist=="HighMed", prune=TRUE)
# phy97.f.c.coral.AD.after.HM <- subset_taxa(phy97.f.c.coral.AD.after.HM, taxa_sums(phy97.f.c.coral.AD.after.HM) > 0, prune=TRUE)
# phy97.f.c.coral.AD.after.HVH <- subset_samples(phy97.f.c.coral.AD.after, Dist=="High" | Dist=="VeryHigh", prune=TRUE)
# phy97.f.c.coral.AD.after.HVH <- subset_taxa(phy97.f.c.coral.AD.after.HVH, taxa_sums(phy97.f.c.coral.AD.after.HVH) > 0, prune=TRUE)

# Calculate number of sequences in phy97.f.c
total_seqs <- sum(taxa_sums(phy97.f.c))

# Cleanup 
rm(a,b,c,i,nam,VeryHigh,VeryLow,phy.f,Low,LowMed,High,HighMed,phy.f.coral)

# Save grouped data as RData file
save(list = ls(all.names = TRUE), file = "data/KI_seqs_f_coral_grouped.RData")
