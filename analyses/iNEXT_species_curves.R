# Clear working directory
rm(list=ls())

# Load necessary packages
if(!require(iNEXT)) { devtools::install_github("JohnsonHsieh/iNEXT"); library(iNEXT) }
library(vegan)
library(gridExtra)
library(metagMisc)

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Color changes
sitecols["27"]<- "#a31b4d"
sitecols["35"]<- "#66a61e"

# Phyloseq to otu helper function 
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Convert phyloseq objects to vegan-friendly format
coral <- vegan_otu(phy97.f.c.coral)
sediment <- vegan_otu(phy97.f.c.sediment)
water <- vegan_otu(phy97.f.c.water)

# Set x limits for each figure
xlim_dist <- c(0,100000)
ylim_dist <- c(0,130)
xlim_coral <- c(0,4000000)
ylim_coral <- c(0,170)
xlim_sed <- c(0,100000)
ylim_sed <- c(0,85)
xlim_wat <- c(0,20000)
ylim_wat <- c(0,85)

# Merge samples by disturbance level
coral_dist <- merge_samples(phy97.f.c.coral,"Dist")
# Extract otu table
coral_dist_otu <- otu_table(coral_dist)
# Rename rownames
rownames(coral_dist_otu) <- c("Medium","Very High")
# Prepare otu table for iNEXT analysis
coral_dist_iNEXT0 <- prepare_inext(coral_dist_otu)
# Run iNEXT (interpolation-extrapolation) analysis
coral_dist_iNEXT <- iNEXT(coral_dist_iNEXT0,q=c(0))
# Make ggplot object from iNEXT analysis
coral_dist_iNEXT_gg <- ggiNEXT(coral_dist_iNEXT,type=1,se=TRUE)+ # Start plotting
  scale_x_continuous(name="Number of Sequences",limits = xlim_coral)+ # Set x axis
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_coral)+ # Set y axis
  theme_classic()+ # Set basic theme
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+ # Set colors
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]])) # Set fill colors
coral_dist_iNEXT_gg 

sed_dist <- merge_samples(phy97.f.c.sediment,"Dist")
sed_dist_otu <- otu_table(sed_dist)
sed_dist_iNEXT0 <- prepare_inext(sed_dist_otu)
sed_dist_iNEXT<- iNEXT(sed_dist_iNEXT0,q=c(0))
sed_dist_iNEXT_gg <- ggiNEXT(sed_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_sed)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_sed)+
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
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))

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

# Create figures
jpeg(filename="figures/coral_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(coral_dist_iNEXT_gg,coral_site_iNEXT_gg,coral_field_season_iNEXT_gg)
dev.off()

jpeg(filename="figures/sediment_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(sed_dist_iNEXT_gg,sed_site_iNEXT_gg,sed_field_season_iNEXT_gg)
dev.off()

jpeg(filename="figures/water_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(wat_dist_iNEXT_gg,wat_site_iNEXT_gg,wat_field_season_iNEXT_gg)
dev.off()

###################### By Coral Species #####################
# Subset coral species
phy97.f.c.coral.Peyd <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Pocillopora_eydouxi")
phy97.f.c.coral.MAeq <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Montipora_foliosa")
phy97.f.c.coral.Plob <- subset_samples(phy97.f.c.coral,sample_data(phy97.f.c.coral)$Coral_Species=="Porites_lobata")

# Set axis limits for plotting
xlim_MAeq <- c(0,1500000)
ylim_MAeq <- c(0,60)
xlim_Peyd <- c(0,1500000)
ylim_Peyd <- c(0,60)
xlim_Plob <- c(0,1000000)
ylim_Plob <- c(0,50)

MAeq_dist <- merge_samples(phy97.f.c.coral.MAeq,"Dist")
MAeq_dist_otu <- otu_table(MAeq_dist)
rownames(MAeq_dist_otu) <- c("Medium","Very High")
MAeq_dist_iNEXT0 <- prepare_inext(MAeq_dist_otu)
MAeq_dist_iNEXT<- iNEXT(MAeq_dist_iNEXT0,q=c(0))
MAeq_dist_iNEXT_gg <- ggiNEXT(MAeq_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_MAeq)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_MAeq)+
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))
MAeq_dist_iNEXT_gg

Peyd_dist <- merge_samples(phy97.f.c.coral.Peyd,"Dist")
Peyd_dist_otu <- otu_table(Peyd_dist)
rownames(Peyd_dist_otu) <- c("Medium","Very High")
Peyd_dist_iNEXT0 <- prepare_inext(Peyd_dist_otu)
Peyd_dist_iNEXT<- iNEXT(Peyd_dist_iNEXT0,q=c(0))
Peyd_dist_iNEXT_gg <- ggiNEXT(Peyd_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_Peyd)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_Peyd)+
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))
Peyd_dist_iNEXT_gg

Plob_dist <- merge_samples(phy97.f.c.coral.Plob,"Dist")
Plob_dist_otu <- otu_table(Plob_dist)
rownames(Plob_dist_otu) <- c("Medium","Very High")
Plob_dist_iNEXT0 <- prepare_inext(Plob_dist_otu)
Plob_dist_iNEXT<- iNEXT(Plob_dist_iNEXT0,q=c(0))
Plob_dist_iNEXT_gg <- ggiNEXT(Plob_dist_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_Plob)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_Plob)+
  theme_classic()+
  scale_color_manual(values=c(sitecols[["8"]],sitecols[["30"]]))+
  scale_fill_manual(values=c(sitecols[["8"]],sitecols[["30"]]))
Plob_dist_iNEXT_gg

###
MAeq_site <- merge_samples(phy97.f.c.coral.MAeq,"site")
MAeq_site_otu <- otu_table(MAeq_site)
MAeq_site_iNEXT0 <- prepare_inext(MAeq_site_otu)
MAeq_site_iNEXT<- iNEXT(MAeq_site_iNEXT0,q=c(0))
MAeq_site_iNEXT_gg <- ggiNEXT(MAeq_site_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_MAeq)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_MAeq)+
  theme_classic()+
  scale_color_manual(values=sitecols)+
  scale_fill_manual(values=sitecols)
MAeq_site_iNEXT_gg

Peyd_site <- merge_samples(phy97.f.c.coral.Peyd,"site")
Peyd_site_otu <- otu_table(Peyd_site)
Peyd_site_iNEXT0 <- prepare_inext(Peyd_site_otu)
Peyd_site_iNEXT<- iNEXT(Peyd_site_iNEXT0,q=c(0))
Peyd_site_iNEXT_gg <- ggiNEXT(Peyd_site_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_Peyd)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_Peyd)+
  theme_classic()+
  scale_color_manual(values=sitecols)+
  scale_fill_manual(values=sitecols)
Peyd_site_iNEXT_gg

Plob_site <- merge_samples(phy97.f.c.coral.Plob,"site")
Plob_site_otu <- otu_table(Plob_site)
Plob_site_iNEXT0 <- prepare_inext(Plob_site_otu)
Plob_site_iNEXT<- iNEXT(Plob_site_iNEXT0,q=c(0))
Plob_site_iNEXT_gg <- ggiNEXT(Plob_site_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_Plob)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_Plob)+
  theme_classic()+
  scale_color_manual(values=sitecols)+
  scale_fill_manual(values=sitecols)
Plob_site_iNEXT_gg

####
MAeq_field_season <- merge_samples(phy97.f.c.coral.MAeq,"field_season")
MAeq_field_season_otu <- otu_table(MAeq_field_season)
MAeq_field_season_iNEXT0 <- prepare_inext(MAeq_field_season_otu)
MAeq_field_season_iNEXT<- iNEXT(MAeq_field_season_iNEXT0,q=c(0))
MAeq_field_season_iNEXT_gg <- ggiNEXT(MAeq_field_season_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_MAeq)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_MAeq)+
  theme_classic()+
  scale_color_manual(values=timecols)+
  scale_fill_manual(values=timecols)
MAeq_field_season_iNEXT_gg

Peyd_field_season <- merge_samples(phy97.f.c.coral.Peyd,"field_season")
Peyd_field_season_otu <- otu_table(Peyd_field_season)
Peyd_field_season_iNEXT0 <- prepare_inext(Peyd_field_season_otu)
Peyd_field_season_iNEXT<- iNEXT(Peyd_field_season_iNEXT0,q=c(0))
Peyd_field_season_iNEXT_gg <- ggiNEXT(Peyd_field_season_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_Peyd)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_Peyd)+
  theme_classic()+
  scale_color_manual(values=timecols)+
  scale_fill_manual(values=timecols)
Peyd_field_season_iNEXT_gg

Plob_field_season <- merge_samples(phy97.f.c.coral.Plob,"field_season")
Plob_field_season_otu <- otu_table(Plob_field_season)
Plob_field_season_iNEXT0 <- prepare_inext(Plob_field_season_otu)
Plob_field_season_iNEXT<- iNEXT(Plob_field_season_iNEXT0,q=c(0))
Plob_field_season_iNEXT_gg <- ggiNEXT(Plob_field_season_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = xlim_Plob)+
  scale_y_continuous(name="OTU Richness",expand=c(0,0),limits = ylim_Plob)+
  theme_classic()+
  scale_color_manual(values=timecols)+
  scale_fill_manual(values=timecols)
Plob_field_season_iNEXT_gg

# Make figures
jpeg(filename="figures/MAeq_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(MAeq_dist_iNEXT_gg,MAeq_site_iNEXT_gg,MAeq_field_season_iNEXT_gg)
dev.off()

jpeg(filename="figures/Peyd_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(Peyd_dist_iNEXT_gg,Peyd_site_iNEXT_gg,Peyd_field_season_iNEXT_gg)
dev.off()

jpeg(filename="figures/Plob_iNEXT_SAC.jpeg",res = 300, width = 9, height= 12,units = "in")
grid.arrange(Plob_dist_iNEXT_gg,Plob_site_iNEXT_gg,Plob_field_season_iNEXT_gg)
dev.off()


######### ALL SAMPLES TOGETHER ########
all <- merge_samples(phy97.f.c,"SampleType")
all_otu <- otu_table(all)
all_iNEXT0 <- prepare_inext(all)
all_iNEXT<- iNEXT(all_iNEXT0,q=c(0))
all_iNEXT_gg <- ggiNEXT(all_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences")+
  scale_y_continuous(name="OTU Richness",expand=c(0,0))+
  theme_classic()+
  theme(legend.position=c(0.85,0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.text = element_text(size=24))+
  scale_color_manual(values=compcols)+
  scale_fill_manual(values=compcols)
all_iNEXT_gg

all_iNEXT_gg2 <- ggiNEXT(all_iNEXT,type=1,se=TRUE)+
  scale_x_continuous(name="Number of Sequences",limits = c(0,100000),breaks = c(seq(0,100000,50000)))+
  scale_y_continuous(name="OTU Richness",expand=c(0,0))+
  theme_classic()+
  theme(legend.position=c(0.8,0.5),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18))+
  scale_color_manual(values=compcols,guide=FALSE)+
  scale_fill_manual(values=compcols, guide=FALSE)+
  guides(fill=FALSE,color=FALSE,shape=FALSE,linetype=FALSE)
all_iNEXT_gg2

# Make inset of zoomed version
vp <- viewport(width = 0.45, height = 0.6, x = 0.44, y = 0.5)
print(all_iNEXT_gg)
print(all_iNEXT_gg2,vp=vp)

# Create figure
jpeg(filename="figures/all_iNEXT_SAC.jpeg",res = 300, width = 9, height= 4,units = "in")
print(all_iNEXT_gg) # Plot the main panel
print(all_iNEXT_gg2,vp=vp) # Plot the inset
dev.off() # Close jpg
