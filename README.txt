main folder

# Files # 
KI_Compartment.Rproj : RProject file for this project
ITS2db_trimmed_derep_data.fasta : ITS2 database used in bioinformatic pipeline

##########

data

# Scripts #

Coralphoto__Metadata/process_coral_metadata.R : append coral photo metadata to mapping file
INPUT Coral Photo Metadata: data/Coralphoto__Metadata/KI_Coralphoto_Metadata_Jan_to_Apr_2017_23March.csv
INPUT Mapping File: data/mapping_file_dada.txt
OUTPUT data/Coralphoto__Metadata/KI_Compartment_metadata.csv
OUTPUT data/Coralphoto__Metadata/KI_Compartment_metadata.csv

export_physeq_objects.R : extract taxa table and otu tables
INPUT analyses/KI_Compartment.RData
OUTPUT Taxa table: data/tax_table.txt
OUTPUT OTU table: data/otu_table.txt

count_water_sediment_samples.R : count number of samples
INPUT data/KI_Compartment_f_coral_grouped.RData
OUTPUT manually extract sample sizes

# Files #
KI_Compartment_f_coral_grouped.RData : all objects created after filtering dataset
KI_Compartment_phyASVfcp.RData : only phyASVfcp after filtering dataset
KI_env_all_KI_Compartment.csv : Environmental parameters from Kiritimati field monitoring program
mapping_file_dada.txt : mapping file for dada2 bioinformatics pipeline

# Folders #
tree/ contains distance matrices, sequences for each genus, and aligne sequences for each genus, as well as the 'uber.tree' which is produced by analyses/dada2/KI_Compartment_filter_samples_dada.R

####################

analyses

# Scripts #

dada2/dada2.R : Dada2 bioinformatic pipeline. Runs the following - 
	- cutadapt to remove primers
	- full dada2 pipeline
INPUT data/Bioinf/sequences/DClaar_2-34439409_seqs_KI_Compartment : original sequences (download from Zenodo to run this script)
INPUT ITS2db_trimmed_derep_dada.fasta : reference ITS2 database
INPUT data/mapping_file_dada.txt : mapping file for sequences
OUTPUT analyses/dada2/KI_Compartment_dada.RData

KI_Compartment_filter_samples_dada.R 
	- filters samples
	- defines site names/disturbance levels
	- make phylogenetic tree
	- transform to proportional abundances
	- count number of sequences
	- subset phyloseqs by compartment/field season
INPUT analyses/dada2/KI_Compartment_dada.RData
OUTPUT data/Bioinf/tree/uber.tre : phylogenetic tree
OUTPUT data/KI_Compartment_f_coral_grouped.RData : all objects created in this script
OUTPUT data/KI_Compartment_phyASVfcp.RData : only phyASVfcp


**environmental_parameters/env_summ
**balance_sensitivity.R
**balance_sensitivity_field_season.R
**count_seqs.R
**count_taxa.R

KI_Compartment_colors.R : creates colors for plots
OUTPUT KI_Compartment_colors.RData 

**KI_Compartment_indicatortaxa.R
**permanova_betadisper_sensitivity.R

# Files #

**environmental_parameters/env_summ.csv
**environmental_parameters/EnvParams.xlsx
KI_Compartment_colors.RData : colors for plots
dada2/KI_Compartment_dada.RData : data after dada2 pipeline, but before additional filtering

# Folders #

**refseqs
**environmental_parameters

#################

# Folders #

**Figure_1
**Figure_2
**Figure_3
**Figure_4
**Figure_5
**Figure_S1
**Figure_S2
**Figure_S3
**Figure_S4
**Figure_S5
**Figure_S6
KI_map

# Files #
**FIgure_2b_coral_water_sed_venn_OTU.pdf
**KI_Compartment_Table_S1.docx
**SupplementaryTable3.xlsx