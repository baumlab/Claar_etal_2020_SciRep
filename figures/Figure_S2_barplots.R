rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(gridExtra)
library(ggplot2)
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


clade_colors <- c("A"="#ffffb3","C"="#8dd3c7","D"="#bebada","F"="#fb8072","G"="#fdb462","I"="#b3de69")


p_MAeq_M_clade <- plot_bar(phy97.f.c.coral.MAeq.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_MAeq_VH_clade <- plot_bar(phy97.f.c.coral.MAeq.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_MAeq_VH_clade

p_Peyd_M_clade <- plot_bar(phy97.f.c.coral.Peyd.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_Peyd_VH_clade <- plot_bar(phy97.f.c.coral.Peyd.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_Peyd_VH_clade

p_Plob_M_clade <- plot_bar(phy97.f.c.coral.Plob.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_Plob_VH_clade <- plot_bar(phy97.f.c.coral.Plob.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_Plob_VH_clade

p_sediment_M_clade <- plot_bar(phy97.f.c.sediment.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_sediment_VH_clade <- plot_bar(phy97.f.c.sediment.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_sediment_VH_clade


p_water_M_clade <- plot_bar(phy97.f.c.water.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Medium Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")+  # Remove ticks 
  coord_flip()
p_water_VH_clade <- plot_bar(phy97.f.c.water.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+coord_flip()
p_water_VH_clade

jpeg(filename = "figures/Fig_S2C_barplot_MAeq.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_MAeq_M_clade,p_MAeq_VH_clade,nrow=2)
dev.off()

jpeg(filename = "figures/Fig_S2D_barplot_Peyd.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_Peyd_M_clade,p_Peyd_VH_clade,nrow=2)
dev.off()

jpeg(filename = "figures/Fig_S2E_barplot_Plob.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_Plob_M_clade,p_Plob_VH_clade,nrow=2)
dev.off()


jpeg(filename = "figures/Fig_S2A_barplot_sediment.jpg",width = 7.5, height = 10,units = "in",res=300)
grid.arrange(p_sediment_M_clade,p_sediment_VH_clade,nrow=2)
dev.off()

jpeg(filename = "figures/Fig_S2B_barplot_water.jpg",width = 7.5, height = 8,units = "in",res=300)
grid.arrange(p_water_M_clade,p_water_VH_clade,nrow=2,heights=c(1.35,1))
dev.off()

# jpeg(filename = "figures/barplot_all.jpg",width = 7.5, height = 10,units = "in",res=300)
# grid.arrange(p_MAeq_M_clade,p_MAeq_VH_clade,p_Peyd_M_clade,p_Peyd_VH_clade,p_Plob_M_clade,p_Plob_VH_clade,p_sediment_M_clade,p_sediment_VH_clade,p_water_M_clade,p_water_VH_clade,nrow=10)
# dev.off()

p_all <- plot_bar(phy97.f.c,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Very High Disturbance")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.direction = "horizontal")+coord_flip()+
  guides(fill=guide_legend(ncol=6))

leg <- get_legend(p_all)
jpeg(filename = "figures/Fig_S2_barplot_clade_legend.jpg",width = 6, height = 1,units = "in",res=300)
as_ggplot(leg)
dev.off()
