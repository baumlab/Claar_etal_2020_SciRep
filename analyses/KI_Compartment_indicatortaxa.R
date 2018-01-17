rm(list=ls())

library(vegan)


load("data/KI_Compartment_f_coral_grouped.RData")

phy97.f.c.sediment <- subset_samples(phy97.f.c.sediment,field_season!="KI2015c")
phy97.f.c <- subset_samples(phy97.f.c,field_season!="KI2015c")

col <- c(KI2014 = "#2c7fb8", KI2015a_Pre = "#7fcdbb", KI2015a_Post = "#253494", KI2015b = "#41b6c4")

# Helper function from http://joey711.github.io/phyloseq-demo/phyloseq-demo.html to convert phyloseq objects to vegan-friendly formatting
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

water <- subset_taxa(physeq = phy97.f.c.water,taxa_sums(phy97.f.c.water) > 50)
sediment <- subset_taxa(physeq = phy97.f.c.sediment,taxa_sums(phy97.f.c.sediment) > 50)


water.veg <- veganotu(water)
sediment.veg <- veganotu(sediment)
all.veg <- veganotu(phy97.f.c)

groups.water <- as.numeric(sample_data(water)$field_season)
indval.water = multipatt(as.data.frame(water.veg), groups.water, control = how(nperm=999))
summary(indval.water)

groups.sed <- as.numeric(sample_data(sediment)$field_season)
indval.sed = multipatt(as.data.frame(sediment.veg), groups.sed, control = how(nperm=999))
summary(indval.sed)

groups.all <- as.numeric(sample_data(phy97.f.c)$SampleType)
indval.all = multipatt(as.data.frame(all.veg), groups.all, control = how(nperm=999))
summary(indval.all)
