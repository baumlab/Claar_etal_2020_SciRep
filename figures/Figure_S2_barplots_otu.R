rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(gridExtra)
library(ggplot2)
library(viridis)
library(phyloseq)
library(ggpubr)

phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

phy97.f.c.coral.MAeq.M.p <- transform_sample_counts(phy97.f.c.coral.MAeq.M, function(x) x/sum(x))
phy97.f.c.coral.MAeq.M.p.f <- subset_taxa(phy97.f.c.coral.MAeq.M.p,taxa_sums(phy97.f.c.coral.MAeq.M.p)>0.001)
phy97.f.c.coral.MAeq.VH.p <- transform_sample_counts(phy97.f.c.coral.MAeq.VH, function(x) x/sum(x))
phy97.f.c.coral.MAeq.VH.p.f <- subset_taxa(phy97.f.c.coral.MAeq.VH.p,taxa_sums(phy97.f.c.coral.MAeq.VH.p)>0.001)
phy97.f.c.coral.Plob.M.p <- transform_sample_counts(phy97.f.c.coral.Plob.M, function(x) x/sum(x))
phy97.f.c.coral.Plob.M.p.f <- subset_taxa(phy97.f.c.coral.Plob.M.p,taxa_sums(phy97.f.c.coral.Plob.M.p)>0.001)
phy97.f.c.coral.Plob.VH.p <- transform_sample_counts(phy97.f.c.coral.Plob.VH, function(x) x/sum(x))
phy97.f.c.coral.Plob.VH.p.f <- subset_taxa(phy97.f.c.coral.Plob.VH.p,taxa_sums(phy97.f.c.coral.Plob.VH.p)>0.001)
phy97.f.c.coral.Peyd.M.p <- transform_sample_counts(phy97.f.c.coral.Peyd.M, function(x) x/sum(x))
phy97.f.c.coral.Peyd.M.p.f <- subset_taxa(phy97.f.c.coral.Peyd.M.p,taxa_sums(phy97.f.c.coral.Peyd.M.p)>0.001)
phy97.f.c.coral.Peyd.VH.p <- transform_sample_counts(phy97.f.c.coral.Peyd.VH, function(x) x/sum(x))
phy97.f.c.coral.Peyd.VH.p.f <- subset_taxa(phy97.f.c.coral.Peyd.VH.p,taxa_sums(phy97.f.c.coral.Peyd.VH.p)>0.001)
phy97.f.c.water.M.p <- transform_sample_counts(phy97.f.c.water.M, function(x) x/sum(x))
phy97.f.c.water.M.p.f <- subset_taxa(phy97.f.c.water.M.p,taxa_sums(phy97.f.c.water.M.p)>0.001)
phy97.f.c.water.VH.p <- transform_sample_counts(phy97.f.c.water.VH, function(x) x/sum(x))
phy97.f.c.water.VH.p.f <- subset_taxa(phy97.f.c.water.VH.p,taxa_sums(phy97.f.c.water.VH.p)>0.001)
phy97.f.c.sediment.M.p <- transform_sample_counts(phy97.f.c.sediment.M, function(x) x/sum(x))
phy97.f.c.sediment.M.p.f <- subset_taxa(phy97.f.c.sediment.M.p,taxa_sums(phy97.f.c.sediment.M.p)>0.001)
phy97.f.c.sediment.VH.p <- transform_sample_counts(phy97.f.c.sediment.VH, function(x) x/sum(x))
phy97.f.c.sediment.VH.p.f <- subset_taxa(phy97.f.c.sediment.VH.p,taxa_sums(phy97.f.c.sediment.VH.p)>0.001)


otus <- rowSums(otu_table(phy97.f.c))
otus <- otus[otus>0]
otu_colors <- c(col=rainbow(351))
denovo_otu <- data.frame(hit=data.frame(tax_table(phy97.f.c))$hit)
rownames(denovo_otu) <- rownames(tax_table(phy97.f.c))
names(otu_colors) <- names(otus)
hit_colors <- c(col=rainbow(399))
names(hit_colors) <- denovo_otu$hit

p_master <- plot_bar(phy97.f.c.p,fill="hit")

p_MAeq_M_otu <- plot_bar(phy97.f.c.coral.MAeq.M.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "bottom")+  # Remove ticks 
  coord_flip()
p_MAeq_VH_otu <- plot_bar(phy97.f.c.coral.MAeq.VH.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_MAeq_VH_otu

p_Peyd_M_otu <- plot_bar(phy97.f.c.coral.Peyd.M.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_Peyd_VH_otu <- plot_bar(phy97.f.c.coral.Peyd.VH.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_Peyd_VH_otu

p_Plob_M_otu <- plot_bar(phy97.f.c.coral.Plob.M.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_Plob_VH_otu <- plot_bar(phy97.f.c.coral.Plob.VH.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_Plob_VH_otu

p_sediment_M_otu <- plot_bar(phy97.f.c.sediment.M.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_sediment_VH_otu <- plot_bar(phy97.f.c.sediment.VH.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_sediment_VH_otu


p_water_M_otu <- plot_bar(phy97.f.c.water.M.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()

p_water_VH_otu <- plot_bar(phy97.f.c.water.VH.p,fill="hit")+
  scale_fill_manual(values=hit_colors,name="hit")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")+coord_flip()
p_water_VH_otu

jpeg(filename = "figures/barplot_MAeq_comb_otu.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_MAeq_M_otu,p_MAeq_VH_otu,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_Peyd_comb_otu.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_Peyd_M_otu,p_Peyd_VH_otu,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_Plob_comb_otu.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_Plob_M_otu,p_Plob_VH_otu,nrow=2)
dev.off()


jpeg(filename = "figures/barplot_sediment_comb_otu.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_sediment_M_otu,p_sediment_VH_otu,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_water_comb_otu.jpg",width = 7.5, height = 8,units = "in",res=300)
grid.arrange(p_water_M_otu,p_water_VH_otu,nrow=2,heights=c(1.35,1))
dev.off()

p_coral_otu2 <- plot_bar(phy97.f.c.coral.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+facet_wrap(~hit)
p_coral_otu2

phy97.f.c.coral.p.01 = filter_taxa(phy97.f.c.coral.p, function(x) mean(x) > 0.01, TRUE)

p_coral_otu3 <- plot_bar(phy97.f.c.coral.p.01,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+facet_wrap(Coral_Species~hit)
p_coral_otu3


phy97.f.c.sediment.p.01 = filter_taxa(phy97.f.c.sediment.p, function(x) mean(x) > 0.01, TRUE)

p_sediment_otu <- plot_bar(phy97.f.c.sediment.p.01,fill="clade")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+facet_wrap(Dist~clade,nrow=2)
p_sediment_otu


gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}
g <- gglegend(p_master)
as_ggplot()

gleg <- get_legend(p_master)
as_ggplot(gleg)



phy97.f.c.s <- phy97.f.c.p
tax_table(phy97.f.c.s)[,8] <- gsub("_.*", "", tax_table(phy97.f.c.s)[,8])

tax_table(phy97.f.c.s) <- tax_table(phy97.f.c.s)[,c(9,8,1)]
phy97.f.c.s.h <- tax_glom(phy97.f.c.s,"hit")

phy97.f.c.s.h
phy97.f.c.s.h <- prune_taxa(taxa_sums(phy97.f.c.s.h)>0,phy97.f.c.s.h)

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


test <- plot_bar(phy97.f.c.s.h.coral.MAeq.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank())
test

plot_bar(phy97.f.c.s.h,fill="hit")

phy97.f.c.s.h.sediment <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$SampleType=="sediment")
phy97.f.c.s.h.sediment <- prune_taxa(taxa_sums(phy97.f.c.s.h.sediment)>0,phy97.f.c.s.h.sediment)
phy97.f.c.s.h.sediment.VH <- subset_samples(phy97.f.c.s.h.sediment,sample_data(phy97.f.c.s.h.sediment)$Dist=="VeryHigh")
phy97.f.c.s.h.sediment.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.sediment.VH)>0,phy97.f.c.s.h.sediment.VH)
phy97.f.c.s.h.sediment.M <- subset_samples(phy97.f.c.s.h.sediment,sample_data(phy97.f.c.s.h.sediment)$Dist=="HighMed")
phy97.f.c.s.h.sediment.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.sediment.M)>0,phy97.f.c.s.h.sediment.M)

p1.sediment <- plot_bar(phy97.f.c.s.h.sediment.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p2.sediment <- plot_bar(phy97.f.c.s.h.sediment.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

grid.arrange(p1.sediment,p2.sediment,nrow=2)

phy97.f.c.s.h.water <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$SampleType=="water")
phy97.f.c.s.h.water <- prune_taxa(taxa_sums(phy97.f.c.s.h.water)>0,phy97.f.c.s.h.water)
phy97.f.c.s.h.water.VH <- subset_samples(phy97.f.c.s.h.water,sample_data(phy97.f.c.s.h.water)$Dist=="VeryHigh")
phy97.f.c.s.h.water.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.water.VH)>0,phy97.f.c.s.h.water.VH)
phy97.f.c.s.h.water.M <- subset_samples(phy97.f.c.s.h.water,sample_data(phy97.f.c.s.h.water)$Dist=="HighMed")
phy97.f.c.s.h.water.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.water.M)>0,phy97.f.c.s.h.water.M)

p1.water <- plot_bar(phy97.f.c.s.h.water.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p2.water <- plot_bar(phy97.f.c.s.h.water.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

grid.arrange(p1.water,p2.water,nrow=2)

phy97.f.c.s.h.coral.Peyd <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.s.h.coral.Peyd <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Peyd)>0,phy97.f.c.s.h.coral.Peyd)
phy97.f.c.s.h.coral.Peyd.VH <- subset_samples(phy97.f.c.s.h.coral.Peyd,sample_data(phy97.f.c.s.h.coral.Peyd)$Dist=="VeryHigh")
phy97.f.c.s.h.coral.Peyd.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Peyd.VH)>0,phy97.f.c.s.h.coral.Peyd.VH)
phy97.f.c.s.h.coral.Peyd.M <- subset_samples(phy97.f.c.s.h.coral.Peyd,sample_data(phy97.f.c.s.h.coral.Peyd)$Dist=="HighMed")
phy97.f.c.s.h.coral.Peyd.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Peyd.M)>0,phy97.f.c.s.h.coral.Peyd.M)

p1.coral.Peyd <- plot_bar(phy97.f.c.s.h.coral.Peyd.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p2.coral.Peyd <- plot_bar(phy97.f.c.s.h.coral.Peyd.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

grid.arrange(p1.coral.Peyd,p2.coral.Peyd,nrow=2)

phy97.f.c.s.h.coral.MAeq <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$Coral_Species=="Montipora_foliosa")
phy97.f.c.s.h.coral.MAeq <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.MAeq)>0,phy97.f.c.s.h.coral.MAeq)
phy97.f.c.s.h.coral.MAeq.VH <- subset_samples(phy97.f.c.s.h.coral.MAeq,sample_data(phy97.f.c.s.h.coral.MAeq)$Dist=="VeryHigh")
phy97.f.c.s.h.coral.MAeq.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.MAeq.VH)>0,phy97.f.c.s.h.coral.MAeq.VH)
phy97.f.c.s.h.coral.MAeq.M <- subset_samples(phy97.f.c.s.h.coral.MAeq,sample_data(phy97.f.c.s.h.coral.MAeq)$Dist=="HighMed")
phy97.f.c.s.h.coral.MAeq.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.MAeq.M)>0,phy97.f.c.s.h.coral.MAeq.M)

p1.coral.MAeq <- plot_bar(phy97.f.c.s.h.coral.MAeq.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p2.coral.MAeq <- plot_bar(phy97.f.c.s.h.coral.MAeq.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

grid.arrange(p1.coral.MAeq,p2.coral.MAeq,nrow=2)

phy97.f.c.s.h.coral.Plob <- subset_samples(phy97.f.c.s.h,sample_data(phy97.f.c.s.h)$Coral_Species=="Porites_lobata")
phy97.f.c.s.h.coral.Plob <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Plob)>0,phy97.f.c.s.h.coral.Plob)
phy97.f.c.s.h.coral.Plob.VH <- subset_samples(phy97.f.c.s.h.coral.Plob,sample_data(phy97.f.c.s.h.coral.Plob)$Dist=="VeryHigh")
phy97.f.c.s.h.coral.Plob.VH <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Plob.VH)>0,phy97.f.c.s.h.coral.Plob.VH)
phy97.f.c.s.h.coral.Plob.M <- subset_samples(phy97.f.c.s.h.coral.Plob,sample_data(phy97.f.c.s.h.coral.Plob)$Dist=="HighMed")
phy97.f.c.s.h.coral.Plob.M <- prune_taxa(taxa_sums(phy97.f.c.s.h.coral.Plob.M)>0,phy97.f.c.s.h.coral.Plob.M)

p1.coral.Plob <- plot_bar(phy97.f.c.s.h.coral.Plob.VH,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p2.coral.Plob <- plot_bar(phy97.f.c.s.h.coral.Plob.M,fill="hit")+
  coord_flip()+
  scale_fill_manual(values=otu_colors2,name="Subclade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

grid.arrange(p1.coral.Plob,p2.coral.Plob,nrow=2)

