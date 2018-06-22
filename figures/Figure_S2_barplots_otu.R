rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(gridExtra)
library(ggplot2)

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

p_MAeq_M_otu <- plot_bar(phy97.f.c.coral.MAeq.M.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_MAeq_VH_otu <- plot_bar(phy97.f.c.coral.MAeq.VH.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_MAeq_VH_otu

p_Peyd_M_otu <- plot_bar(phy97.f.c.coral.Peyd.M.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_Peyd_VH_otu <- plot_bar(phy97.f.c.coral.Peyd.VH.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_Peyd_VH_otu

p_Plob_M_otu <- plot_bar(phy97.f.c.coral.Plob.M.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_Plob_VH_otu <- plot_bar(phy97.f.c.coral.Plob.VH.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_Plob_VH_otu

p_sediment_M_otu <- plot_bar(phy97.f.c.sediment.M.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_sediment_VH_otu <- plot_bar(phy97.f.c.sediment.VH.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_sediment_VH_otu


p_water_M_otu <- plot_bar(phy97.f.c.water.M.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_water_VH_otu <- plot_bar(phy97.f.c.water.VH.p,fill="otu")+
  # scale_fill_manual(values=otu_colors,name="otu")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_water_VH_otu

jpeg(filename = "figures/barplot_MAeq_comb_otu.jpg",width = 7.5, height = 6,units = "in",res=300)
grid.arrange(p_MAeq_M_otu,p_MAeq_VH_otu,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_Peyd_comb_otu.jpg",width = 7.5, height = 6,units = "in",res=300)
grid.arrange(p_Peyd_M_otu,p_Peyd_VH_otu,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_Plob_comb_otu.jpg",width = 7.5, height = 6,units = "in",res=300)
grid.arrange(p_Plob_M_otu,p_Plob_VH_otu,nrow=2)
dev.off()


jpeg(filename = "figures/barplot_sediment_comb_otu.jpg",width = 7.5, height = 6,units = "in",res=300)
grid.arrange(p_sediment_M_otu,p_sediment_VH_otu,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_water_comb_otu.jpg",width = 7.5, height = 6,units = "in",res=300)
grid.arrange(p_water_M_otu,p_water_VH_otu,nrow=2)
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
