rm(list=ls())
if(!require(iNEXT)) { devtools::install_github("JohnsonHsieh/iNEXT"); library(iNEXT) }
library(vegan)
library(gridExtra)
library(metagMisc)
library(phyloseq)
library(ggplot2)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

coral <- vegan_otu(phy97.f.c.coral)
sediment <- vegan_otu(phy97.f.c.sediment)
water <- vegan_otu(phy97.f.c.water)

xlim_dist <- c(0,100000)
ylim_dist <- c(0,130)
xlim_coral <- c(0,4000000)
ylim_coral <- c(0,170)
xlim_sed <- c(0,100000)
ylim_sed <- c(0,85)
xlim_wat <- c(0,20000)
ylim_wat <- c(0,85)

coral_dist <- merge_samples(phy97.f.c.coral,"Dist")
coral_dist_otu <- otu_table(coral_dist)
rownames(coral_dist_otu) <- c("Medium","Very High")
coral_dist_iNEXT0 <- prepare_inext(coral_dist_otu)
coral_dist_iNEXT<- iNEXT(coral_dist_iNEXT0,q=c(0))
coral_dist_iNEXT_gg <- ggiNEXT(coral_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_coral)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_coral)+
  # scale_x_continuous(name="Number of Sequences")+
  # scale_y_continuous(name="OTU Richness",expand=c(0,0))+
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))
coral_dist_iNEXT_gg

sed_dist <- merge_samples(phy97.f.c.sediment,"Dist")
sed_dist_otu <- otu_table(sed_dist)
sed_dist_iNEXT0 <- prepare_inext(sed_dist_otu)
sed_dist_iNEXT<- iNEXT(sed_dist_iNEXT0,q=c(0))
sed_dist_iNEXT_gg <- ggiNEXT(sed_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_sed)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_sed)+
  # scale_x_continuous(name="Number of Sequences")+
  # scale_y_continuous(name="OTU Richness",expand=c(0,0))+
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))

wat_dist <- merge_samples(phy97.f.c.water,"Dist")
wat_dist_otu <- otu_table(wat_dist)
wat_dist_iNEXT0 <- prepare_inext(wat_dist_otu)
wat_dist_iNEXT<- iNEXT(wat_dist_iNEXT0,q=c(0))
wat_dist_iNEXT_gg <- ggiNEXT(wat_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_wat)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_wat)+
  # scale_x_continuous(name="Number of Sequences")+
  # scale_y_continuous(name="OTU Richness",expand=c(0,0))+
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))

sitecols["27"]<- "#a31b4d"
sitecols["35"]<- "#66a61e"

coral_site <- merge_samples(phy97.f.c.coral,"site")
coral_site_otu <- otu_table(coral_site)
coral_site_iNEXT0 <- prepare_inext(coral_site_otu)
coral_site_iNEXT<- iNEXT(coral_site_iNEXT0,q=c(0))
coral_site_iNEXT_gg <- ggiNEXT(coral_site_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_coral)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_coral)+
  theme_classic()+
  scale_color_manual(values=c(sitecols))+
  scale_fill_manual(values=c(sitecols))

sed_site <- merge_samples(phy97.f.c.sediment,"site")
sed_site_otu <- otu_table(sed_site)
sed_site_iNEXT0 <- prepare_inext(sed_site_otu)
sed_site_iNEXT<- iNEXT(sed_site_iNEXT0,q=c(0))
sed_site_iNEXT_gg <- ggiNEXT(sed_site_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_sed)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_sed)+
  theme_classic()+
  scale_color_manual(values=c(sitecols))+
  scale_fill_manual(values=c(sitecols))

wat_site <- merge_samples(phy97.f.c.water,"site")
wat_site_otu <- otu_table(wat_site)
wat_site_iNEXT0 <- prepare_inext(wat_site_otu)
wat_site_iNEXT<- iNEXT(wat_site_iNEXT0,q=c(0))
wat_site_iNEXT_gg <- ggiNEXT(wat_site_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_wat)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_wat)+
  theme_classic()+
  scale_color_manual(values=c(sitecols))+
  scale_fill_manual(values=c(sitecols))

coral_field_season <- merge_samples(phy97.f.c.coral,"field_season")
coral_field_season_otu <- otu_table(coral_field_season)
coral_field_season_iNEXT0 <- prepare_inext(coral_field_season_otu)
coral_field_season_iNEXT<- iNEXT(coral_field_season_iNEXT0,q=c(0))
coral_field_season_iNEXT_gg <- ggiNEXT(coral_field_season_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_coral)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_coral)+
  theme_classic()+
  scale_color_manual(values=c(timecols))+
  scale_fill_manual(values=c(timecols))

sed_field_season <- merge_samples(phy97.f.c.sediment,"field_season")
sed_field_season_otu <- otu_table(sed_field_season)
sed_field_season_iNEXT0 <- prepare_inext(sed_field_season_otu)
sed_field_season_iNEXT<- iNEXT(sed_field_season_iNEXT0,q=c(0))
sed_field_season_iNEXT_gg <- ggiNEXT(sed_field_season_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_sed)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_sed)+
  theme_classic()+
  scale_color_manual(values=c(timecols))+
  scale_fill_manual(values=c(timecols))

wat_field_season <- merge_samples(phy97.f.c.water,"field_season")
wat_field_season_otu <- otu_table(wat_field_season)
wat_field_season_iNEXT0 <- prepare_inext(wat_field_season_otu)
wat_field_season_iNEXT<- iNEXT(wat_field_season_iNEXT0,q=c(0))
wat_field_season_iNEXT_gg <- ggiNEXT(wat_field_season_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_wat)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_wat)+
  theme_classic()+
  scale_color_manual(values=c(timecols))+
  scale_fill_manual(values=c(timecols))


jpeg(filename="figures/Fig_S6_coral_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(coral_dist_iNEXT_gg,coral_site_iNEXT_gg,coral_field_season_iNEXT_gg)
dev.off()

jpeg(filename="figures/Fig_S7_sediment_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(sed_dist_iNEXT_gg,sed_site_iNEXT_gg,sed_field_season_iNEXT_gg)
dev.off()

jpeg(filename="figures/Fig_S8_water_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(wat_dist_iNEXT_gg,wat_site_iNEXT_gg,wat_field_season_iNEXT_gg)
dev.off()

pdf(file="figures/Fig_S6_coral_iNEXT_SAC.pdf",width = 9, height = 12,useDingbats = FALSE)
grid.arrange(coral_dist_iNEXT_gg,coral_site_iNEXT_gg,coral_field_season_iNEXT_gg)
dev.off()

pdf(file="figures/Fig_S7_sediment_iNEXT_SAC.pdf",width = 9, height = 12,useDingbats = FALSE)
grid.arrange(sed_dist_iNEXT_gg,sed_site_iNEXT_gg,sed_field_season_iNEXT_gg)
dev.off()

pdf(file="figures/Fig_S8_water_iNEXT_SAC.pdf", width = 9, height = 12,useDingbats = FALSE)
grid.arrange(wat_dist_iNEXT_gg,wat_site_iNEXT_gg,wat_field_season_iNEXT_gg)
dev.off()

