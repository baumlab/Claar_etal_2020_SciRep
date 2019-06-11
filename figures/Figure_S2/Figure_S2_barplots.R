# Clear working environment
rm(list=ls())

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Load necessary packages
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(phyloseq)

# Subset coral species
phyASV.f.c.coral.Peyd <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phyASV.f.c.coral.MAeq <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Montipora_foliosa")
phyASV.f.c.coral.Plob <- subset_samples(phyASV.f.c.coral,sample_data(phyASV.f.c.coral)$Coral_Species=="Porites_lobata")

# Transform sample counts to percentage
phyASV.f.c.coral.MAeq.M <- subset_samples(phyASV.f.c.coral.MAeq,sample_data(phyASV.f.c.coral.MAeq)$Dist =="Medium")
phyASV.f.c.coral.MAeq.M.p <- transform_sample_counts(phyASV.f.c.coral.MAeq.M, function(x) x/sum(x))
phyASV.f.c.coral.MAeq.M.p.f <- subset_taxa(phyASV.f.c.coral.MAeq.M.p,taxa_sums(phyASV.f.c.coral.MAeq.M.p)>0.001)

phyASV.f.c.coral.MAeq.VH <- subset_samples(phyASV.f.c.coral.MAeq,sample_data(phyASV.f.c.coral.MAeq)$Dist =="VeryHigh")
phyASV.f.c.coral.MAeq.VH.p <- transform_sample_counts(phyASV.f.c.coral.MAeq.VH, function(x) x/sum(x))
phyASV.f.c.coral.MAeq.VH.p.f <- subset_taxa(phyASV.f.c.coral.MAeq.VH.p,taxa_sums(phyASV.f.c.coral.MAeq.VH.p)>0.001)

phyASV.f.c.coral.Plob.M <- subset_samples(phyASV.f.c.coral.Plob,sample_data(phyASV.f.c.coral.Plob)$Dist =="Medium")
phyASV.f.c.coral.Plob.M.p <- transform_sample_counts(phyASV.f.c.coral.Plob.M, function(x) x/sum(x))
phyASV.f.c.coral.Plob.M.p.f <- subset_taxa(phyASV.f.c.coral.Plob.M.p,taxa_sums(phyASV.f.c.coral.Plob.M.p)>0.001)

phyASV.f.c.coral.Plob.VH <- subset_samples(phyASV.f.c.coral.Plob,sample_data(phyASV.f.c.coral.Plob)$Dist =="VeryHigh")
phyASV.f.c.coral.Plob.VH.p <- transform_sample_counts(phyASV.f.c.coral.Plob.VH, function(x) x/sum(x))
phyASV.f.c.coral.Plob.VH.p.f <- subset_taxa(phyASV.f.c.coral.Plob.VH.p,taxa_sums(phyASV.f.c.coral.Plob.VH.p)>0.001)

phyASV.f.c.coral.Peyd.M <- subset_samples(phyASV.f.c.coral.Peyd,sample_data(phyASV.f.c.coral.Peyd)$Dist =="Medium")
phyASV.f.c.coral.Peyd.M.p <- transform_sample_counts(phyASV.f.c.coral.Peyd.M, function(x) x/sum(x))
phyASV.f.c.coral.Peyd.M.p.f <- subset_taxa(phyASV.f.c.coral.Peyd.M.p,taxa_sums(phyASV.f.c.coral.Peyd.M.p)>0.001)

phyASV.f.c.coral.Peyd.VH <- subset_samples(phyASV.f.c.coral.Peyd,sample_data(phyASV.f.c.coral.Peyd)$Dist =="VeryHigh")
phyASV.f.c.coral.Peyd.VH.p <- transform_sample_counts(phyASV.f.c.coral.Peyd.VH, function(x) x/sum(x))
phyASV.f.c.coral.Peyd.VH.p.f <- subset_taxa(phyASV.f.c.coral.Peyd.VH.p,taxa_sums(phyASV.f.c.coral.Peyd.VH.p)>0.001)

phyASV.f.c.water.M <- subset_samples(phyASV.f.c.water,sample_data(phyASV.f.c.water)$Dist =="Medium")
phyASV.f.c.water.M.p <- transform_sample_counts(phyASV.f.c.water.M, function(x) x/sum(x))
phyASV.f.c.water.M.p.f <- subset_taxa(phyASV.f.c.water.M.p,taxa_sums(phyASV.f.c.water.M.p)>0.001)

phyASV.f.c.water.VH <- subset_samples(phyASV.f.c.water,sample_data(phyASV.f.c.water)$Dist =="VeryHigh")
phyASV.f.c.water.VH.p <- transform_sample_counts(phyASV.f.c.water.VH, function(x) x/sum(x))
phyASV.f.c.water.VH.p.f <- subset_taxa(phyASV.f.c.water.VH.p,taxa_sums(phyASV.f.c.water.VH.p)>0.001)

phyASV.f.c.sediment.M <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$Dist =="Medium")
phyASV.f.c.sediment.M <- subset_samples(phyASV.f.c.sediment.M,sample_sums(phyASV.f.c.sediment.M)>0)
phyASV.f.c.sediment.M.p <- transform_sample_counts(phyASV.f.c.sediment.M, function(x) x/sum(x))
phyASV.f.c.sediment.M.p.f <- subset_taxa(phyASV.f.c.sediment.M.p,taxa_sums(phyASV.f.c.sediment.M.p)>0.001)

phyASV.f.c.sediment.VH <- subset_samples(phyASV.f.c.sediment,sample_data(phyASV.f.c.sediment)$Dist =="VeryHigh")
phyASV.f.c.sediment.VH.p <- transform_sample_counts(phyASV.f.c.sediment.VH, function(x) x/sum(x))
phyASV.f.c.sediment.VH.p.f <- subset_taxa(phyASV.f.c.sediment.VH.p,taxa_sums(phyASV.f.c.sediment.VH.p)>0.001)

# Create barplot for Montipora at Medium disturbance sites
p_MAeq_M_genus <- plot_bar(phyASV.f.c.coral.MAeq.M.p,fill="Genus")+ # Start plotting
  scale_fill_manual(values=genus_colors,name="Genus")+ # Set fill colors
  scale_y_continuous(expand = c(0,0))+ # remove extra space between data and axes
  scale_x_discrete(name="Medium Disturbance")+ # Label axis
  theme(axis.text.x = element_blank(), # Remove x axis text
        axis.ticks = element_blank(), # Remove x axis ticks
        axis.title.x = element_blank(), # Remove x axis title
        axis.text.y = element_blank(), # Remove y axis text
        legend.position = "none")+  # Suppress plotting the legend
  coord_flip() # Flip from vertical to horizontal
# Create barplot for Montipora at Very High disturbance sites
p_MAeq_VH_genus <- plot_bar(phyASV.f.c.coral.MAeq.VH.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()

# Pocillopora
p_Peyd_M_genus <- plot_bar(phyASV.f.c.coral.Peyd.M.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+ 
  coord_flip()
p_Peyd_VH_genus <- plot_bar(phyASV.f.c.coral.Peyd.VH.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()

# Porites
p_Plob_M_genus <- plot_bar(phyASV.f.c.coral.Plob.M.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  
  coord_flip()
p_Plob_VH_genus <- plot_bar(phyASV.f.c.coral.Plob.VH.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()

# Sediment
p_sediment_M_genus <- plot_bar(phyASV.f.c.sediment.M.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  
  coord_flip()
p_sediment_VH_genus <- plot_bar(phyASV.f.c.sediment.VH.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()

# Water
p_water_M_genus <- plot_bar(phyASV.f.c.water.M.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  
  coord_flip()
p_water_VH_genus <- plot_bar(phyASV.f.c.water.VH.p,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()

# Make jpegs
jpeg(filename = "figures/Figure_S2/Fig_S2C_barplot_MAeq.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_MAeq_M_genus,p_MAeq_VH_genus,nrow=2)
dev.off()

jpeg(filename = "figures/Figure_S2/Fig_S2D_barplot_Peyd.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_Peyd_M_genus,p_Peyd_VH_genus,nrow=2)
dev.off()

jpeg(filename = "figures/Figure_S2/Fig_S2E_barplot_Plob.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_Plob_M_genus,p_Plob_VH_genus,nrow=2)
dev.off()

jpeg(filename = "figures/Figure_S2/Fig_S2A_barplot_sediment.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_sediment_M_genus,p_sediment_VH_genus,nrow=2)
dev.off()

jpeg(filename = "figures/Figure_S2/Fig_S2B_barplot_water.jpg",width = 7.5, height = 8,units = "in",res=300)
grid.arrange(p_water_M_genus,p_water_VH_genus,nrow=2,heights=c(1.35,1))
dev.off()



# Make a ggplot object for plotting the legend
p_all <- plot_bar(phyASV.f.c,fill="Genus")+
  scale_fill_manual(values=genus_colors,name="Genus")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.direction = "horizontal")+coord_flip()+
  guides(fill=guide_legend(ncol=6))

leg <- get_legend(p_all) # Extract only the legend for plotting
# Create figure
jpeg(filename = "figures/Figure_S2/Fig_S2_barplot_genus_legend.jpg",width = 6, height = 1,units = "in",res=300)
as_ggplot(leg) # Plot the legend
dev.off() # Close jpg

pdf(file = "figures/Figure_S2/Fig_S2_barplot_genus_legend.pdf",width = 6, height = 1)
as_ggplot(leg) # Plot the legend
dev.off() # Close jpg


jpeg(filename = "figures/Figure_S2/Fig_S2_barplot_all.jpg",width = 8.5, height = 11,units = "in",res=300)
grid.arrange(p_Peyd_M_genus+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(name="P. grandis")+ggtitle("Medium"), 
             p_Peyd_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle("Very High"),
             p_MAeq_M_genus+scale_x_discrete(name="M. aequituberculata"),
             p_MAeq_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
             p_Plob_M_genus+scale_x_discrete(name="P. lobata"), 
             p_Plob_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
             p_sediment_M_genus+scale_x_discrete(name="Sediment"), 
             p_sediment_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
             p_water_M_genus+theme(axis.title.x = element_text("Abundance"))+scale_x_discrete(name="Water"), 
             p_water_VH_genus+theme(axis.title.y = element_blank()),
             nrow=5,ncol=2)
dev.off()

pdf(file = "figures/Figure_S2/Fig_S2_barplot_all.pdf",width = 8.5, height = 11)
grid.arrange(p_Peyd_M_genus+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(name="P. grandis")+ggtitle("Medium"), 
             p_Peyd_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle("Very High"),
             p_MAeq_M_genus+scale_x_discrete(name="M. aequituberculata"),
             p_MAeq_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
             p_Plob_M_genus+scale_x_discrete(name="P. lobata"), 
             p_Plob_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
             p_sediment_M_genus+scale_x_discrete(name="Sediment"), 
             p_sediment_VH_genus+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
             p_water_M_genus+theme(axis.title.x = element_text("Abundance"))+scale_x_discrete(name="Water"), 
             p_water_VH_genus+theme(axis.title.y = element_blank()),
             nrow=5,ncol=2)
dev.off()