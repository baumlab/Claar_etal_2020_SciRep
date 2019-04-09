# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Load necessary packages
library(gridExtra)
library(ggplot2)
library(viridis)
library(phyloseq)
library(ggpubr)

# Create new phyloseq object for this figure
phy97.f.c.s <- phy97.f.c.p
# Remove accession number information, only keep subclade identifier in "hit"
tax_table(phy97.f.c.s)[,8] <- gsub("_.*", "", tax_table(phy97.f.c.s)[,8])

sample_data(phy97.f.c.s)$Dist <- gsub("27","VeryHigh",sample_data(phy97.f.c.s)$Dist)

# Reorder tax table so that tax glom actually works
tax_table(phy97.f.c.s) <- tax_table(phy97.f.c.s)[,c(9,8,1)]
# Agglomerate by subclade
phy97.f.c.s.h <- tax_glom(phy97.f.c.s,"hit")
# Renive any taxa that are empty (sum=0)
phy97.f.c.s.h <- prune_taxa(taxa_sums(phy97.f.c.s.h)>0,phy97.f.c.s.h)

# Make a color palette for plotting
otus2 <- rowSums(otu_table(phy97.f.c.s.h))
otus2 <- otus2[otus2>0]
otu_colors2 <- c(col=rainbow(length(otus2)))
denovo_otu2 <- data.frame(hit=data.frame(tax_table(phy97.f.c.s.h))$hit)
rownames(denovo_otu2) <- rownames(tax_table(phy97.f.c.s.h))
names(otu_colors2) <- denovo_otu2$hit

# Set specific colors so that things are actually distinguishable
otu_colors2["G100"] <- "burlywood4"
otu_colors2["C1.v1a"] <- "cyan3"
otu_colors2["C1d"] <- "cornflowerblue"
otu_colors2["I2"] <- "#bc7a10"
otu_colors2["I4"] <- "#d18408"
otu_colors2["C15"] <- "aquamarine4"
otu_colors2["C15.4"] <- "aquamarine2"
otu_colors2["C1226"] <- "aquamarine3"
otu_colors2["C1"] <- "#AEFF00"
otu_colors2["C3"] <- "#00FF93"
otu_colors2["C31a"] <- "#4DE800"
otu_colors2["C32.1"] <- "#9AFFAE"
otu_colors2["C42"] <- "#2A7B2D"
otu_colors2["G3.3"] <- "#51260d"
otu_colors2["C31"] <- "#9FC9A1"




# Subset samples for each figure
# Sediment
phy97.f.c.s.h.sediment <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$SampleType=="sediment")
phy97.f.c.s.h.sediment <- prune_taxa(taxa_sums(phy97.f.c.s.h.sediment)>0,phy97.f.c.s.h.sediment)
phy97.f.c.s.h.sediment.VH <- subset_samples(phy97.f.c.s.h.sediment,sample_data(phy97.f.c.s.h.sediment)$Dist=="VeryHigh")
phy97.f.c.s.h.sediment.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.sediment.VH)>0,phy97.f.c.s.h.sediment.VH)
phy97.f.c.s.h.sediment.M <- subset_samples(phy97.f.c.s.h.sediment,sample_data(phy97.f.c.s.h.sediment)$Dist=="Medium")
phy97.f.c.s.h.sediment.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.sediment.M)>0,phy97.f.c.s.h.sediment.M)

# Water
phy97.f.c.s.h.water <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$SampleType=="water")
phy97.f.c.s.h.water <- prune_taxa(taxa_sums(phy97.f.c.s.h.water)>0,phy97.f.c.s.h.water)
phy97.f.c.s.h.water.VH <- subset_samples(phy97.f.c.s.h.water,sample_data(phy97.f.c.s.h.water)$Dist=="VeryHigh")
phy97.f.c.s.h.water.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.water.VH)>0,phy97.f.c.s.h.water.VH)
phy97.f.c.s.h.water.M <- subset_samples(phy97.f.c.s.h.water,sample_data(phy97.f.c.s.h.water)$Dist=="Medium")
phy97.f.c.s.h.water.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.water.M)>0,phy97.f.c.s.h.water.M)

# Pocillopora
phy97.f.c.s.h.coral.Peyd <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.s.h.coral.Peyd <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Peyd)>0,phy97.f.c.s.h.coral.Peyd)
phy97.f.c.s.h.coral.Peyd.VH <- subset_samples(phy97.f.c.s.h.coral.Peyd,sample_data(phy97.f.c.s.h.coral.Peyd)$Dist=="VeryHigh")
phy97.f.c.s.h.coral.Peyd.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Peyd.VH)>0,phy97.f.c.s.h.coral.Peyd.VH)
phy97.f.c.s.h.coral.Peyd.M <- subset_samples(phy97.f.c.s.h.coral.Peyd,sample_data(phy97.f.c.s.h.coral.Peyd)$Dist=="Medium")
phy97.f.c.s.h.coral.Peyd.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Peyd.M)>0,phy97.f.c.s.h.coral.Peyd.M)

# Montipora
phy97.f.c.s.h.coral.MAeq <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$Coral_Species=="Montipora_foliosa")
phy97.f.c.s.h.coral.MAeq <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.MAeq)>0,phy97.f.c.s.h.coral.MAeq)
phy97.f.c.s.h.coral.MAeq.VH <- subset_samples(phy97.f.c.s.h.coral.MAeq,sample_data(phy97.f.c.s.h.coral.MAeq)$Dist=="VeryHigh")
phy97.f.c.s.h.coral.MAeq.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.MAeq.VH)>0,phy97.f.c.s.h.coral.MAeq.VH)
phy97.f.c.s.h.coral.MAeq.M <- subset_samples(phy97.f.c.s.h.coral.MAeq,sample_data(phy97.f.c.s.h.coral.MAeq)$Dist=="Medium")
phy97.f.c.s.h.coral.MAeq.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.MAeq.M)>0,phy97.f.c.s.h.coral.MAeq.M)

# Porites
phy97.f.c.s.h.coral.Plob <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$Coral_Species=="Porites_lobata")
phy97.f.c.s.h.coral.Plob <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Plob)>0,phy97.f.c.s.h.coral.Plob)
phy97.f.c.s.h.coral.Plob.VH <- subset_samples(phy97.f.c.s.h.coral.Plob,sample_data(phy97.f.c.s.h.coral.Plob)$Dist=="VeryHigh")
phy97.f.c.s.h.coral.Plob.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Plob.VH)>0,phy97.f.c.s.h.coral.Plob.VH)
phy97.f.c.s.h.coral.Plob.M <- subset_samples(phy97.f.c.s.h.coral.Plob,sample_data(phy97.f.c.s.h.coral.Plob)$Dist=="Medium")
phy97.f.c.s.h.coral.Plob.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Plob.M)>0,phy97.f.c.s.h.coral.Plob.M)

# Bar plots
# Sediment - Very High
p1.sediment <- plot_bar(phy97.f.c.s.h.sediment.VH,fill="hit")+ # start plotting
  coord_flip()+ # flip horizontal
  scale_fill_manual(values=otu_colors2,name="Subclade")+ # set color fill values
  scale_y_continuous(expand = c(0, 0))+ # remove extra space between data and axes
  scale_x_discrete(name="Very High Disturbance")+ # rename x axis
  theme(axis.text.x = element_blank(), # remove x axis text
        axis.ticks = element_blank(), # remove x axis ticks
        axis.text.y = element_blank(), # remove y axis text
        legend.position = "none") # suppress plotting of legend 
# Sediment - Medium
p2.sediment <- plot_bar(phy97.f.c.s.h.sediment.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.title.x = element_blank(), # remove x axis title
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# Water - Very High
p1.water <- plot_bar(phy97.f.c.s.h.water.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")
# Water - Medium
p2.water <- plot_bar(phy97.f.c.s.h.water.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# Pocillopora - Very High
p1.coral.Peyd <- plot_bar(phy97.f.c.s.h.coral.Peyd.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")
# Pocillopora - Medium
p2.coral.Peyd <- plot_bar(phy97.f.c.s.h.coral.Peyd.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# Montipora - Very High
p1.coral.MAeq <- plot_bar(phy97.f.c.s.h.coral.MAeq.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")
# Montipora - Medium
p2.coral.MAeq <- plot_bar(phy97.f.c.s.h.coral.MAeq.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# Porites - Very High
p1.coral.Plob <- plot_bar(phy97.f.c.s.h.coral.Plob.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")
# Porites - Medium
p2.coral.Plob <- plot_bar(phy97.f.c.s.h.coral.Plob.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# Plot Very High and Medium together 
grid.arrange(p1.sediment,p2.sediment,nrow=2)
grid.arrange(p1.water,p2.water,nrow=2)
grid.arrange(p1.coral.Peyd,p2.coral.Peyd,nrow=2)
grid.arrange(p1.coral.MAeq,p2.coral.MAeq,nrow=2)
grid.arrange(p1.coral.Plob,p2.coral.Plob,nrow=2)

dev.off()
# Make jpegs
# Sediment
jpeg(filename = "figures/Fig_S3B_barplot_sediment_subclade.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p2.sediment,p1.sediment,nrow=2)
dev.off()
# Water
jpeg(filename = "figures/Fig_S3C_barplot_water_subclade.jpg",width = 7.5, height = 8,units = "in",res=300)
grid.arrange(p2.water,p1.water,nrow=2,heights=c(1.35,1))
dev.off()
# Pocillopora
jpeg(filename = "figures/Fig_S3E_barplot_Peyd_subclade.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p2.coral.Peyd,p1.coral.Peyd,nrow=2)
dev.off()
# Montipora
jpeg(filename = "figures/Fig_S3D_barplot_MAeq_subclade.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p2.coral.MAeq,p1.coral.MAeq,nrow=2)
dev.off()
# Porites
jpeg(filename = "figures/Fig_S3F_barplot_Plob_subclade.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p2.coral.Plob,p1.coral.Plob,nrow=2)
dev.off()

# # For testing colors  
# test <- plot_bar(phy97.f.c.s.h.coral.MAeq.M,fill="hit")+
#   coord_flip()+
#   scale_fill_manual(values=otu_colors2,name="Subclade")+
#   scale_y_continuous(expand = c(0, 0))+
#   scale_x_discrete(name="Very High Disturbance")+
#   theme(axis.text.x = element_blank(), # Remove x axis tick labels
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank())
# test

# Make ggplot with all for legend creation
p0 <- plot_bar(phy97.f.c.s.h,fill="hit")+  
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  theme(legend.position = "bottom")

# Get the legend
leg <- get_legend(p0)
# Make a figure with just the legend
jpeg(filename = "figures/Fig_S3A_barplot_subclade_legend.jpg",width = 5.5, height = 5.5,units = "in",res=300)
as_ggplot(leg) # Plot the legend
dev.off() # Close jpg

### Calculate number of dominant subclades for each compartment
otu1 <- otu_table(phy97.f.c.s.h.sediment)
tax1 <- tax_table(phy97.f.c.s.h.sediment)
otu_tax1 <- cbind(otu1,tax1)
max_sub <- apply(otu1, 2, function(x) which(x == max(x)))
sort(unique(max_sub))
nsub1 <- length(unique(max_sub))

otu2 <- otu_table(phy97.f.c.s.h.water)
tax2 <- tax_table(phy97.f.c.s.h.water)
otu_tax2 <- cbind(otu2,tax2)
max_sub2 <- apply(otu2, 2, function(x) which(x == max(x)))
sort(unique(max_sub2))
nsub2 <- length(unique(max_sub2))

otu3 <- otu_table(phy97.f.c.s.h.coral.Plob)
tax3 <- tax_table(phy97.f.c.s.h.coral.Plob)
otu_tax3 <- cbind(otu3,tax3)
max_sub3 <- apply(otu3, 2, function(x) which(x == max(x)))
sort(unique(max_sub3))
nsub3 <- length(unique(max_sub3))

### Calculate total number of  subclades for each compartment
ntaxa(phy97.f.c.s.h.sediment)
ntaxa(phy97.f.c.s.h.water)
ntaxa(phy97.f.c.s.h.coral.MAeq)
ntaxa(phy97.f.c.s.h.coral.Peyd)
ntaxa(phy97.f.c.s.h.coral.Plob)

unique(data.frame(tax_table(phy97.f.c.s.h.sediment))$clade)
unique(data.frame(tax_table(phy97.f.c.s.h.water))$clade)
unique(data.frame(tax_table(phy97.f.c.s.h.coral.MAeq))$clade)
unique(data.frame(tax_table(phy97.f.c.s.h.coral.Peyd))$clade)
unique(data.frame(tax_table(phy97.f.c.s.h.coral.Plob))$clade)
