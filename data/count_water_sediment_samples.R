rm(list=ls())
load("data/KI_Compartment_f_coral_grouped.RData")

phyASV.f.c.water

w <- data.frame(sample_data(phyASV.f.c.water))
s <- data.frame(sample_data(phyASV.f.c.sediment))

nrow(w[w$site=="8"&w$field_season=="KI2014",])
nrow(w[w$site=="8"&w$field_season=="KI2015a_Post",])
nrow(w[w$site=="8"&w$field_season=="KI2015b",])

nrow(w[w$site=="27"&w$field_season=="KI2014",])
nrow(w[w$site=="27"&w$field_season=="KI2015a_Post",])
nrow(w[w$site=="27"&w$field_season=="KI2015b",])

nrow(w[w$site=="30"&w$field_season=="KI2014",])
nrow(w[w$site=="30"&w$field_season=="KI2015a_Post",])
nrow(w[w$site=="30"&w$field_season=="KI2015b",])

nrow(w[w$site=="35"&w$field_season=="KI2014",])
nrow(w[w$site=="35"&w$field_season=="KI2015a_Post",])
nrow(w[w$site=="35"&w$field_season=="KI2015b",])

# Sediment
nrow(s[s$site=="8"&s$field_season=="KI2014",])
nrow(s[s$site=="8"&s$field_season=="KI2015a_Post",])
nrow(s[s$site=="8"&s$field_season=="KI2015b",])

nrow(s[s$site=="27"&s$field_season=="KI2014",])
nrow(s[s$site=="27"&s$field_season=="KI2015a_Post",])
nrow(s[s$site=="27"&s$field_season=="KI2015b",])

nrow(s[s$site=="30"&s$field_season=="KI2014",])
nrow(s[s$site=="30"&s$field_season=="KI2015a_Post",])
nrow(s[s$site=="30"&s$field_season=="KI2015b",])

nrow(s[s$site=="35"&s$field_season=="KI2014",])
nrow(s[s$site=="35"&s$field_season=="KI2015a_Post",])
nrow(s[s$site=="35"&s$field_season=="KI2015b",])
