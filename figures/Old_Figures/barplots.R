rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

library(gridExtra)
library(ggplot2)

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
  scale_x_discrete(name="Local Disturbance = Medium")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        legend.position = "none")   # Remove ticks 

p_MAeq_VH_clade <- plot_bar(phy97.f.c.coral.MAeq.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Very High")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())   

p_Plob_M_clade <- plot_bar(phy97.f.c.coral.Plob.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Medium")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        legend.position = "none")   # Remove ticks 

p_Plob_VH_clade <- plot_bar(phy97.f.c.coral.Plob.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Very High")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())   

p_Peyd_M_clade <- plot_bar(phy97.f.c.coral.Peyd.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Medium")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        legend.position = "none")   # Remove ticks 

p_Peyd_VH_clade <- plot_bar(phy97.f.c.coral.Peyd.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Very High")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())   

p_water_M_clade <- plot_bar(phy97.f.c.water.M.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Medium")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        legend.position = "none")   # Remove ticks 

p_water_VH_clade <- plot_bar(phy97.f.c.water.VH.p,fill="clade")+
  scale_fill_manual(values=clade_colors,name="Clade")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(name="Local Disturbance = Very High")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())   

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

jpeg(filename = "figures/barplot_sediment_comb.jpg",width = 7.5, height = 6,units = "in",res=300)
grid.arrange(p_sediment_M_clade,p_sediment_VH_clade,nrow=2)
dev.off()

jpeg(filename = "figures/barplot_sediment.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_sediment_M_clade,p_sediment_VH_clade,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Sediment",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_water.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_water_M_clade,p_water_VH_clade,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Water",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_MAeq.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_MAeq_M_clade,p_MAeq_VH_clade,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Montipora aequituberculata",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_Peyd.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_Peyd_M_clade,p_Peyd_VH_clade,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Pocillopora grandis",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_Plob.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_Plob_M_clade,p_Plob_VH_clade,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Porites lobata",gp=gpar(fontsize=20,font=3)))
dev.off()











gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



p_MAeq_M_hit <- plot_bar(phy97.f.c.coral.MAeq.M.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")


p_MAeq_VH_hit <- plot_bar(phy97.f.c.coral.MAeq.VH.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")

p_Plob_M_hit <- plot_bar(phy97.f.c.coral.Plob.M.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")


p_Plob_VH_hit <- plot_bar(phy97.f.c.coral.Plob.VH.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")

p_Peyd_M_hit <- plot_bar(phy97.f.c.coral.Peyd.M.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")


p_Peyd_VH_hit <- plot_bar(phy97.f.c.coral.Peyd.VH.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")

p_water_M_hit <- plot_bar(phy97.f.c.water.M.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")


p_water_VH_hit <- plot_bar(phy97.f.c.water.VH.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")

p_sediment_M_hit <- plot_bar(phy97.f.c.sediment.M.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")


p_sediment_VH_hit <- plot_bar(phy97.f.c.sediment.VH.p,fill="hit")+
  theme(axis.text.x = element_blank(), # Remove x axis tick labels
        legend.position="bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_blank(),
        axis.ticks = element_blank())+   # Remove ticks 
  scale_y_continuous(expand = c(0, 0))+guides(fill=guide_legend(ncol=4))+
  scale_x_discrete(name="Sample")+
  scale_fill_manual(values = hit.col)
p_sediment_VH_hit


hit.col <- c(gg_color_hue(113))
# hit.col <- hue_pal()(113)
hit.col <- data.frame(color=c(hit.col[50:113],hit.col[1:49]))
hit.col <- data.frame(color=hit.col)
hit.col$hit <- unique(data.frame(tax_table(phy97.f.c))$hit)

hit.col <- setNames(as.character(hit.col$color), hit.col$hit)
# Have to figure out how to do this so that it will use the same colors for the same hit on every plot. I think I'm super close, but it isn't quite working right now....

jpeg(filename = "figures/barplot_sediment_hit.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_sediment_M_hit,p_sediment_VH_hit,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Sediment",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_water_hit.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_water_M_hit,p_water_VH_hit,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Water",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_MAeq_hit.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_MAeq_M_hit,p_MAeq_VH_hit,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Montipora aequituberculata",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_Peyd_hit.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_Peyd_M_hit,p_Peyd_VH_hit,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Pocillopora grandis",gp=gpar(fontsize=20,font=3)))
dev.off()

jpeg(filename = "figures/barplot_Plob_hit.jpg",width = 7.5, height = 4,units = "in",res=300)
grid.arrange(p_Plob_M_hit,p_Plob_VH_hit,
             ncol=2,widths=c(1,1.075),
             top = textGrob("Porites lobata",gp=gpar(fontsize=20,font=3)))
dev.off()

