library(dichromat)
library(maptools)
library(scales)
library(RColorBrewer)
library(rgdal)
library(ggplot2)

### site data
sites<-read.csv('figures/KI_map/KI_sites_compartment.csv')

###village data
villages<-read.csv("figures/KI_map/KI_villagesDCC_2015update.csv", header = TRUE) # you get an error but it works

villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.66 + 0.005, 1.69 + 0.005, 1.71 + 0.005))
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.5522183, -157.5522183, -157.5522183)
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"


## set palette for fishing pressure
fishing.cols<-c("#8c510a","#c7eae5")
fishing.cols<-as.data.frame(fishing.cols)
fishing.cols$f.pressure<-levels(sites$f.pressure)
sites$col<-fishing.cols$fishing.cols[match(sites$f.pressure, fishing.cols$f.pressure)]

jpeg(file="figures/Figure1.jpg",width = 7.6, height = 7.2,units="in",res=300)
source("figures/KI_map/KI_base_B&W.R")
# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, 
        add=TRUE,inches=0.3, bg=alpha("black", 0.5))
## legend for village size
text(-157.435 - 0.008, 1.66, "Village with 1500 people", cex=0.6)   
text(-157.435 - 0.008, 1.688, "Village with 1000 people", cex=0.6)   
text(-157.4385 - 0.008, 1.709, "Village with 500 people", cex=0.6)   
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4) 
with(sites, text(lon, lat, label=site, cex=0.5))
legend(-157.57, 1.81,legend=levels(sites$f.pressure), pt.bg=c("#c7eae5","#8c510a"), pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.59, 1.78, "Human\ndisturbance", srt=90, cex=0.6)
text(-157.585, 1.685, "Village", srt=90, cex=0.6)
segments(-157.569, 1.73,-157.569, 1.83)
segments(-157.569, 1.65,-157.569, 1.72)

dev.off()

