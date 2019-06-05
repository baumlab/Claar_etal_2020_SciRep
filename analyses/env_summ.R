rm(list=ls())

library(rlang)
library(skimr)

env <- read.csv("data/KI_env_all_KI_Compartment.csv")

env <- env[env$field_season=="2014"| env$field_season=="2015b",]

names(env)

skim(env)

env_skimmed <- env %>%
  dplyr::group_by(site) %>%
  skim()

env_skimmed 

str(env_skimmed)

env_skimmed2 <- env %>%
  dplyr::group_by(site) %>%
  skim_to_list()
env_skimmed3 <- env %>%
  dplyr::group_by(fpressure) %>%
  skim_to_list()

env_skimmed4 <- env %>%
  dplyr::group_by(fpressure) %>%
  skim_to_wide()

env_skimmed2[["numeric"]]

env_skimmed2[["numeric"]] %>%  
  dplyr::select(site, variable, mean, sd) %>%
  dplyr::filter(variable == "nitrate_plus_nitrite_uM_mean")
env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "nitrate_plus_nitrite_uM_mean")

env_skimmed2[["numeric"]] %>%  
  dplyr::select(site, variable, mean, sd) %>%
  dplyr::filter(variable == "visibility_m")
env_skim_vis <- env_skimmed2[["numeric"]] %>%  
  dplyr::select(site, variable, mean, sd) %>%
  dplyr::filter(variable == "visibility_m")
env_skim_vis2 <- env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "visibility_m")
write.csv(x=env_skim_vis2,file="analyses/env_skim_vis.csv",row.names = FALSE)


env_skimmed2[["numeric"]] %>%  
  dplyr::select(site, variable, mean, sd) %>%
  dplyr::filter(variable == "ysi_salinity_mean_1m")
env_skim_sal <- env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "ysi_salinity_mean_1m")
write.csv(x=env_skim_sal,file="analyses/env_skim_sal.csv",row.names = FALSE)


env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "turb_mean_aq")

env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "phosphate_uM_mean")

env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "prod_mean_f")

env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "prod_mean_aq")

env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "npp_mean_sat")

env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "npp_max_sat")

env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd) %>%
  dplyr::filter(variable == "ysi_DO_mean_1m")

env_skim_all <- env_skimmed3[["numeric"]] %>%  
  dplyr::select(fpressure, variable, mean, sd)
write.csv(x=env_skim_all,file="analyses/env_skim_all.csv",row.names = FALSE)


vis_anova <- aov(visibility_m~fpressure+field_season,data=env)
summary(vis_anova) 
nit_anova <- aov(nitrate_plus_nitrite_uM_mean~fpressure+field_season,data=env)
summary(nit_anova) 
phos_anova <- aov(phosphate_uM_mean~fpressure+field_season,data=env)
summary(phos_anova) 
prod_aq_anova <- aov(prod_mean_aq~fpressure+field_season,data=env)
summary(prod_aq_anova) 
prod_f_anova <- aov(prod_mean_f~fpressure+field_season,data=env)
summary(prod_f_anova) 
turb_f_anova <- aov(turb_mean_aq~fpressure+field_season,data=env)
summary(turb_f_anova) 
DO_anova <- aov(ysi_DO_mean_1m~fpressure,data=env)
summary(DO_anova) 
DO_mgL_anova <- aov(ysi_DO_mgL_mean_1m~fpressure,data=env)
summary(DO_mgL_anova) 
pH_anova <- aov(ysi_pH_mean_1m~fpressure,data=env)
summary(pH_anova) 
npp_sat_anova <- aov(npp_mean_sat~fpressure+field_season,data=env)
summary(npp_sat_anova) 
npp_min_sat_anova <- aov(npp_min_sat~fpressure,data=env)
summary(npp_min_sat_anova) 
npp_max_sat_anova <- aov(npp_max_sat~fpressure,data=env)
summary(npp_max_sat_anova) 
wave_sat_anova <- aov(wave_mean_sat~fpressure,data=env)
summary(wave_sat_anova) 
wave_wind_anova <- aov(wave_wind_fetch_sat~fpressure,data=env)
summary(wave_wind_anova) 


TukeyHSD(vis_anova) # Medium is higher than very high
TukeyHSD(nit_anova) # Medium is higher than very high
TukeyHSD(phos_anova) # Medium is higher than very high
TukeyHSD(prod_aq_anova) # Very High is higher than Medium
TukeyHSD(DO_anova) # Medium is higher than Very High
TukeyHSD(npp_sat_anova) # Medium is higher than Very High
TukeyHSD(npp_min_sat_anova) # Medium is higher than Very High
TukeyHSD(npp_max_sat_anova) # Very High is higher than Medium
TukeyHSD(wave_sat_anova) # Medium is higher than Very High
TukeyHSD(wave_wind_anova) # Very High is higher than Medium

env$NtoP <- env$nitrate_plus_nitrite_uM_mean/env$phosphate_uM_mean

NtoP_anova <- aov(NtoP~fpressure,data=env)
summary(NtoP_anova) 


## Set which metrics to include in PCA
include <- c("visibility_m", "nitrate_plus_nitrite_uM_mean", 
             "phosphate_uM_mean", "prod_mean_aq", "turb_mean_aq", 
             "ysi_DO_mean_1m","NtoP")

env.summ <- env %>%
  group_by(site) %>%
  summarise_at(vars(include), funs(mean(., na.rm=TRUE)))

env.summ2 <- env %>%
  group_by(fpressure) %>%
  summarise_at(vars(include), funs(mean(., na.rm=TRUE)))

include2 <- c("visibility_m", "nitrate_plus_nitrite_uM_mean", 
             "phosphate_uM_mean", "prod_mean_aq", "prod_mean_f",
             "turb_mean_aq", "ysi_temp_mean_1m", "ysi_salinity_mean_1m",
             "ysi_DO_mean_1m","ysi_pH_mean_1m", 
             "npp_mean_sat","npp_max_sat","wave_mean_sat")
env.summ3 <- env %>%
  group_by(fpressure) %>%
  summarise_at(vars(include2), funs(mean(., na.rm=TRUE), sd(., na.rm=TRUE)))
write.csv(x=env.summ3,file="analyses/env_summ.csv",row.names = FALSE)


## Subset temperature metrics for PCA
env_incl <- as.data.frame(env.summ[, include])

## Run PCA
pca.env <- prcomp(na.omit(env_incl), center=TRUE, scale. = TRUE)
par(mfrow=c(1,1))
screeplot(pca.env)

## Calculate loading values of original metrics
aload <- abs(pca.env$rotation)
sweep(aload, 2, colSums(aload), "/")

## Extract temperature PC values for each site
env.PCs <- pca.env$x

## Visualize PCA as biplot
biplot(pca.env)

summary(pca.env)


load("analyses/KI_Compartment_colors.RData")


# Calculate ranges of each column
env2 <- rbind(apply(as.data.frame(env.summ), 2, function(x) rev(range(x))),
               as.data.frame(env.summ))
# Reorder columns
env2 <- env2[, c(2:7)]
# Plot radar chart
rad <- radarchart(
  env2,
  #custom polygon
  plwd=4, plty=1,
  #custom colors
  pcol=sitecols[c(1,2,2,1)],
  pfcol=scales::alpha(sitecols[c(1,2,2,1)], 0.25),
  #custom the grid
  cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
  #custom labels
  vlcex=0.8 
)

legend("topright", legend=env.summ$site, lty=1, lwd=2, col=sitecols[c(1,2,2,1)], bty="n")
