library(dichromat)
library(maptools)
library(scales)
library(RColorBrewer)
library(rgdal)
library(ggplot2)

rm(list=ls())

### site data
sites<-read.csv('figures/KI_map/KI_sites_compartment.csv')
sites$site.simple <- c("M1","VH2","VH1","M2")

###village data
villages<-read.csv("figures/KI_map/KI_villagesDCC_2015update.csv", header = TRUE)

villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.66 - 0.008, 1.69 - 0.008, 1.71 - 0.008))
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.5522183 + 0.01, -157.5522183 + 0.01, -157.5522183 + 0.01)
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"


## set palette for fishing pressure
fishing.cols<-c("#c7eae5","#8c510a")
fishing.cols<-as.data.frame(fishing.cols)
fishing.cols$f.pressure<-levels(sites$f.pressure)
sites$col<-fishing.cols$fishing.cols[match(sites$f.pressure, fishing.cols$f.pressure)]

# jpeg(filename = "figures/KI_map.jpg",width = 7, height = 7.2, units="in", res=300)
jpeg(filename = "figures/Figure_1_KI_map.jpg",width = 7, height = 7.2, units="in", res=300)
source("figures/KI_map/KI_base_B&W.R")
# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("black", 0.4))
## legend for village size
text(-157.478 + 0.01, 1.648, "1500 people", cex=0.69)   
text(-157.478 + 0.01, 1.6735, "1000 people", cex=0.69) 
text(-157.482 + 0.01, 1.699, "500 people", cex=0.69) 
text(-157.588 + 0.004, 1.67, "Village", srt=90, cex=0.6)
segments(-157.563, 1.64,-157.563, 1.704)  
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4)
with(sites, text(lon, lat, label=site.simple, cex=0.3))
legend(-157.56, 1.8,legend=levels(sites$f.pressure), pt.bg=c("#c7eae5","#8c510a"), pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.588, 1.77, "Human\nDisturbance", srt=90, cex=0.6)
segments(-157.563, 1.72,-157.563, 1.82)
# text(-157.3, 1.88, "Bay of\nWrecks", cex = 0.5)
# text(-157.53, 1.82, "Vaskess\nBay", cex = 0.5)

dev.off()

# jpeg(filename = "figures/KI_map.jpg",width = 7, height = 7.2, units="in", res=300)
pdf(file = "figures/Figure_1_KI_map.pdf",width = 7, height = 7.2, useDingbats = FALSE)
source("figures/KI_map/KI_base_B&W.R")
# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("black", 0.4))
## legend for village size
text(-157.478 + 0.01, 1.648, "1500 people", cex=0.69)   
text(-157.478 + 0.01, 1.6735, "1000 people", cex=0.69) 
text(-157.482 + 0.01, 1.699, "500 people", cex=0.69) 
text(-157.588 + 0.004, 1.67, "Village", srt=90, cex=0.6)
segments(-157.563, 1.64,-157.563, 1.704)  
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4)
with(sites, text(lon, lat, label=site.simple, cex=0.3))
legend(-157.56, 1.8,legend=levels(sites$f.pressure), pt.bg=c("#c7eae5","#8c510a"), pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.588, 1.77, "Human\nDisturbance", srt=90, cex=0.6)
segments(-157.563, 1.72,-157.563, 1.82)

dev.off()

