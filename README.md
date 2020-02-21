
## Data and R code accompanying:  
  
**Chronic disturbance modulates symbiont (Symbiodiniaceae) beta diversity on a coral reef**  
  
Authors: Danielle C. Claar, Kristina L. Tietjen, Kieran D. Cox, Ruth D. Gates, Julia K. Baum

_______  


# Structure of files and folders within this repo 

## main folder

### Files
KI\_Compartment.Rproj : RProject file for this project  
ITS2db\_trimmed\_derep\_data.fasta : ITS2 database used in bioinformatic pipeline

## data

### Scripts

**Coralphoto\_\_Metadata/process\_coral\_metadata.R** : append coral photo metadata to mapping file  
INPUT Coral Photo Metadata: data/Coralphoto\_\_Metadata/KI\_Coralphoto\_Metadata\_Jan\_to\_Apr\_2017\_23March.csv  
INPUT Mapping File: data/mapping\_file\_dada.txt  
OUTPUT  data/Coralphoto\_\_Metadata/KI\_Compartment\_metadata.csv  
OUTPUT  data/Coralphoto\_\_Metadata/KI\_Compartment\_metadata.csv  

**export\_physeq\_objects.R** : extract taxa table and otu tables  
INPUT analyses/KI\_Compartment.RData  
OUTPUT Taxa table: data/tax\_table.txt  
OUTPUT OTU table: data/otu\_table.txt  

**count_water_sediment_samples.R** : count number of samples  
INPUT data/KI\_Compartment\_f\_coral_grouped.RData  
OUTPUT manually extract sample sizes  

### Files 
**KI\_Compartment\_f\_coral\_grouped.RData** : all objects created after filtering dataset  
**KI\_Compartment\_phyASVfcp.RData** : only phyASVfcp after filtering dataset  
**KI\_env\_all\_KI\_Compartment.csv** : Environmental parameters from Kiritimati field monitoring program  
**mapping\_file\_dada.txt** : mapping file for dada2 bioinformatics pipeline  

### Folders
**tree/** contains distance matrices, sequences for each genus, and aligne sequences for each genus, as well as the 'uber.tree' which is produced by analyses/dada2/KI\_Compartment\_filter\_samples\_dada.R  

## analyses

### Scripts

**dada2/dada2.R** : Dada2 bioinformatic pipeline. Runs the following -  
- cutadapt to remove primers  
- full dada2 pipeline  
INPUT data/Bioinf/sequences/DClaar\_2-34439409\_seqs\_KI\_Compartment : original sequences (download from Zenodo to run this script)  
INPUT ITS2db_trimmed_derep_dada.fasta : reference ITS2 database  
INPUT data/mapping_file_dada.txt : mapping file for sequences  
OUTPUT analyses/dada2/KI_Compartment_dada.RData  
  
**dada2/KI\_Compartment\_filter\_samples_dada.R**   
- filters samples  
- defines site names/disturbance levels  
- make phylogenetic tree  
- transform to proportional abundances  
- count number of sequences  
- subset phyloseqs by compartment/field season  
INPUT analyses/dada2/KI\_Compartment\_dada.RData  
OUTPUT data/Bioinf/tree/uber.tre : phylogenetic tree  
OUTPUT data/KI\_Compartment\_f\_coral_grouped.RData : all objects created in this script  
OUTPUT data/KI\_Compartment\_phyASVfcp.RData : only phyASVfcp  

**environmental\_parameters/env\_summ.R** : summarize environmental parameters from field surveys and satellite data  
INPUT data/KI\_env\_all\_KI\_Compartment.csv  
OUTPUT analyses/env\_summ.csv  

**count\_seqs.R** : count number of sequences (total, and by compartment)  

**count\_taxa.R** : count number of taxa by coral species  

**KI\_Compartment\_colors.R** : creates colors for plots  
OUTPUT KI\_Compartment\_colors.RData   

**balance\_sensitivity.R** : sensitivity analysis for Disturbance, to see if PERMANOVA results remain significant if dataset is pared down to a balanced design (see methods for justification)  

**balance\_sensitivity\_field\_season.R** : sensitivity analysis for Field Season, to see if PERMANOVA results remain significant if dataset is pared down to a balanced design (see methods for justification)

**permanova\_betadisper\_sensitivity.R** : calculate number of samples to pare down to for "balance\_sensitivity*.R" scripts

### Files

**environmental\_parameters/env\_summ.csv** : summary of environmental parameters  
**KI\_Compartment\_colors.RData** : colors for plots
**dada2/KI\_Compartment\_dada.RData** : data after dada2 pipeline, but before additional filtering

### Folders

**refseqs/** : dominant sequences in .fasta format for blasting. Also includes refseq\_blast\_IDs.xlsx, with blast information

## figures ##

### Folders  
Each figure folder contains the script to produce the figure, the final figure, and any intermediate files  

Figure\_1  
Figure\_2  
Figure\_3  
Figure\_4  
Figure\_5  
Figure\_S1  
Figure\_S2  
Figure\_S3  
Figure\_S4  
Figure\_S5  
Figure\_S6  

### Files  
KI\_Compartment\_Table\_S1.docx  
SupplementaryTable3.xlsx  