#!/bin/bash

find . -type f -exec gunzip {} + # Unzip these files

# They are going to get really big! And this may take quite some time!

python data/Bioinf/make_configs.py # Make the Illumina-utils config files for the Bokulich QC method

python data/Bioinf/boku_qc.py # Run Bokulich QC method on all applicable files
chmod u+x data/Bioinf/boku_qc_KI_Compartment.sh


sh # Give permissions so .sh file will run
source ~/virtual-envs/illumina-utils-v2.0.0/bin/activate
data/Bioinf/boku_qc_KI_Compartment.sh # Run .sh file that was just created 

python data/Bioinf/make_merge_configs.py # Make the Illumina-utils config files for merging pairs

python data/Bioinf/iu_merge_pairs.py # Merge pairs using Illumina-utils
chmod u+x data/Bioinf/iu-merge-pairs_KI_Compartment.sh # Give permissions so .sh file will run
data/Bioinf/iu-merge-pairs_KI_Compartment.sh # Run .sh file that was just created 
# find . -maxdepth 1 -name # To check how far along the merging is (because it takes a long time to do hundreds of samples)

python data/Bioinf/iu_filter_merged_reads.py # Filter merged reads (MAX-MISMATCH=3) using Illumina-utils
chmod u+x data/Bioinf/iu-filter-merged-reads_KI_Compartment.sh # Give permissions so .sh file will run
data/Bioinf/iu-filter-merged-reads_KI_Compartment.sh # Run .sh file that was just created 

python data/Bioinf/rename.py # Rename merged files for downstream processing

