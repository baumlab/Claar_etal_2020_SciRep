# Clear working directory
rm(list=ls())

# Load necessary packages
library(rlang)
library(skimr)
library(dplyr)

# Load data
env <- read.csv("data/KI_env_all_KI_Compartment.csv")

# Limit to field seasons 2014 and 2015b
env <- env[env$field_season=="2014"| env$field_season=="2015b",]
names(env)
# skim(env)
# env_skimmed3 <- env %>%
#   dplyr::group_by(fpressure) %>%
#   skim_to_list()

# env_skimmed4 <- env %>%
#   dplyr::group_by(fpressure) %>%
#   skim_to_wide()

include2 <- c("visibility_m", "ysi_temp_mean_1m", "ysi_salinity_mean_1m",
             "ysi_DO_mean_1m","ysi_pH_mean_1m", 
             "npp_mean_sat","npp_max_sat","wave_mean_sat")
env.summ3 <- env %>%
  group_by(fpressure) %>%
  summarise_at(vars(include2), funs(mean(., na.rm=TRUE), sd(., na.rm=TRUE)))
write.csv(x=env.summ3,file="analyses/env_summ.csv",row.names = FALSE)