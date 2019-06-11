# KI Compartment Colors

# Clear the working directory
rm(list=ls())

# Colors for each field season
timecols <- c(KI2014 = "#2c7fb8", KI2015a_Pre = "#7fcdbb", KI2015a_Post = "#253494", KI2015b = "#41b6c4")

# Colors for each site
sitecols <- c("8" = "#1b9e77","30"="#d95f02","35" ="#7570b3", "27"="#66a61e")

# Colors for each compartment
compcols <- c(coral = "#bf812d", sediment = "#d8b365", water = "#5ab4ac")

# Colors for each coral species
speccols <- c("M. aequituberculata" = "#ec7014", "P. lobata" = "#cc7282", "P. eydouxi" = "#fe9929")

# Colors for each Symbiodinium clade
clade_colors <- c("A"="#ffffb3","C"="#8dd3c7","D"="#bebada","F"="#fb8072","G"="#fdb462","I"="#b3de69")

# Colors for each Symbiodiniaceae genus
genus_colors <- c("g__Symbiodinium"="#ffffb3","g__Breviolum"="gray","g__Cladocopium"="#8dd3c7","g__Durusdinium"="#bebada",
                  "g__Fugacium"="#fb8072","g__Gerakladium"="#fdb462","g__cladeI"="#b3de69")

# Save 
save.image(file="analyses/KI_Compartment_colors.RData")


