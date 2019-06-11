# Load necessary libraries and programs
library(dada2)
library(phyloseq)
# Make sure you have access to cutadapt on your machine
# See https://cutadapt.readthedocs.io/en/stable/installation.html for installation information
cutadapt <- "/usr/bin/cutadapt" # Where is cutadapt located on your machine
system2(cutadapt, args = "--version") # Check that cutadapt works

# https://benjjneb.github.io/dada2/ITS_workflow.html

# Create file lists
path <- "data/Bioinf/sequences/DClaar_2-34439409_seqs_KI_Compartment" # Path to where sequences are stored
fnFs <- sort(list.files(path, pattern="_R1.fastq"))	
fnRs <- sort(list.files(path, pattern="_R2.fastq"))	

fnFs <- file.path(path, fnFs) # Append the filepath to the filenames - forward reads
fnRs <- file.path(path, fnRs) # Append the filepath to the filenames - reverse reads
fnFs
fnRs

# Identify primers used for this project
FWD <- "GTGAATTGCAGAACTCCGTG" # Note that these are only the ITS2 primers (do not include the Illumina adapters, which are already removed)
REV <- "CCTCCGCTTACTTATATGCTT" # Note that these are only the ITS2 primers (do not include the Illumina adapters, which are already removed)

# Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Pre-filter the sequences just to remove those with Ns, but perform no other filtering for now
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, # dada does not allow Ns, so must be zero
              multithread = TRUE, 
              compress = FALSE) # Turn off compression (default is TRUE), cutadapt doesn't work on compressed files

# Count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
# Check to make sure that something shows up for both forward and reverse - this can be a good "spell check", the first time I did it I didn't get any in forward, due to a typo.

# Use cutadapt to remove primers
path.cut <- file.path(path, "cutadapt") # Append the filepath to the filenames - cutadapt
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs)) # Append the filepath to the filenames - forward reads for cutadapt
fnRs.cut <- file.path(path.cut, basename(fnRs)) # Append the filepath to the filenames - reverse reads for cutadapt

FWD.RC <- dada2:::rc(FWD) # Create object that is reverse complement of forward reads
REV.RC <- dada2:::rc(REV) # Create object that is reverse complement of reverse reads
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Sanity check to make sure that cutadapt worked and all primers are gone
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
### There were still some left after the first time, FWD.ReverseReads (Forward) = 6 and REV.ForardReads (Forward) = 21, 
# so I tried to run cutadapt again on the output, but it didn't change anything

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Look at quality profiles of trimmed reads
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

# Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# ref_seqs <- readDNAStringSet(filepath = "ITS2db_fromSymPortal2.fasta")
ref_seqs <- readDNAStringSet(filepath = "ITS2db_trimmed_derep_dada.fasta")
min(width(ref_seqs))
max(width(ref_seqs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, # Standard filtering parameter
                     maxEE = c(2, 2), # Sets the maximum number of “expected errors” allowed in a read
                     trimRight = 100, # The number of nucleotides to remove from the end of each read - chosen based on declining quality scors in plot
                     minLen = 150, # Enforce a minLen here, to get rid of spurious short sequences; was 50 before
                     rm.phix = TRUE, # Remove any phiX sequences
                     compress = TRUE, 
                     multithread = TRUE, # on windows, set multithread = FALSE
                     verbose = TRUE) 
head(out)

# Learn error rates
# Please ignore all the “Not all sequences were the same length.” messages in the next couple sections. We know they aren’t, and it’s OK!
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates as a sanity check
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference - apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      trimOverhang=TRUE, verbose=TRUE)

# Construct amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab.nochim)))

# Track reads through the pipeline - inspect the # of reads that made it through each step in the pipeline to verify everything worked as expected.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
sym.ref <- "ITS2db_trimmed_derep_dada.fasta" 
sym.ref2 <- "ITS2db_fromSymPortal2.fasta" 

# Try with in-house reference database down to Genus
taxa <- assignTaxonomy(seqtab.nochim, sym.ref, multithread = TRUE, tryRC = TRUE, minBoot = 80, verbose = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
tp <- data.frame(taxa.print)
tp$Genus
unique(tp$Genus)

# # Try with SymPortal database
# taxa2 <- assignTaxonomy(seqtab.nochim, sym.ref2, multithread = TRUE, tryRC = TRUE, minBoot = 80, verbose = TRUE)
# taxa.print2 <- taxa2  # Removing sequence rownames for display only
# rownames(taxa.print2) <- NULL
# head(taxa.print2)
# tp2 <- data.frame(taxa.print2)
# tp2$Genus
# unique(tp2$Genus)

# Import to phyloseq
samdf <- read.table("data/mapping_file_dada.txt",header = TRUE) # Read in sample data
rownames(samdf) <- samdf$SampleID
head(samdf)

# We now construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# It is more convenient to use short names for our ASVs (e.g. ASV21) rather than the full DNA sequence 
# when working with some of the tables and visualizations from phyloseq, but we want to keep the full DNA 
# sequences for other purposes like merging with other datasets or indexing into reference databases like the 
# Earth Microbiome Project. For that reason we’ll store the DNA sequences of our ASVs in the refseq slot of the 
# phyloseq object, and then rename our taxa to a short string. That way, the short new taxa names will appear 
# in tables and plots, and we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps).

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

save(ps, file="analyses/dada2/KI_Compartment_dada.RData")

