
## Script for creating KI map for platy ms fall 2017
setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")

library(dichromat)
library(maptools)
library(scales)
library(RColorBrewer)
library(rgdal)
library(ggplot2)

### site data
sites<-read.csv('ki_map_files/ki_sites_platy.csv')

###village data
villages<-read.csv("ki_map_files/KI_villagesDCC_2015update.csv", header = TRUE) # you get an error but it works

villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.66, 1.69, 1.71))
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.5522183, -157.5522183, -157.5522183)
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"



### reordering levels for colouring plot
levels(sites$f.pressure)

sites$levels<-as.numeric(factor(sites$f.pressure, levels(factor(sites$f.pressure))[c(4,1,3,2,5)]))
sites$f.pressure<-factor(sites$f.pressure, levels(factor(sites$f.pressure))[c(4,1,3,2,5)])

## set palette for fishing pressure

fishing.cols<-c("#8c510a","#d8b365","#c7eae5","#5ab4ac","#01665e")
fishing.cols<-as.data.frame(fishing.cols)
fishing.cols$f.pressure<-levels(sites$f.pressure)
sites$col<-fishing.cols$fishing.cols[match(sites$f.pressure, fishing.cols$f.pressure)]


### Sites sampled for platy paper
setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")
#pdf(file="KI_map_sites_platyms.pdf")
#tiff(file="KI_map_sites_platyms.tiff",width = 7, height = 7,units="in",res=300)
#jpeg(file="KI_map_sites_platyms.jpeg",width = 7, height = 7,units="in",res=300)
source("ki_map_files/KI_base_B&W.R")
points(sites$lon, sites$lat, bg=alpha(sites$col,0.7), pch=21, cex=3.0)
#with(sites_platy, text(lon, lat, label=site, cex=0.7))
legend(-157.23, 2.08,legend=levels(sites$f.pressure), pt.bg=c("#8c510a","#d8b365","#c7eae5","#5ab4ac","#01665e"), pch=21, bty="n", pt.cex=2.0, cex=1)
text(-157.24, 2.035, "Human disturbance", srt=90, cex=1)
segments(-157.23, 2.004,-157.23, 2.068)
text(-157.3, 1.9, "Bay\n of\n Wrecks", cex=1)
text(-157.52, 1.82, "Vaskess\n Bay", cex=1)

dev.off()

############### ok going to try to put villages on the same map as the sites ####################
setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")
#tiff(file="KI_map_platysites_villages.tiff",width = 7, height = 7,units="in",res=300)
#jpeg(file="KI_map_platysites_villages.jpeg",width = 7, height = 7,units="in",res=300)
#pdf(file="KI_map_platysites_villages.pdf")
source("ki_map_files/KI_base_B&W.R")

# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("red", 0.4))
## legend for village size
text(-157.465, 1.66, "Village with 1500 people", cex=1)   
text(-157.465, 1.688, "Village with 1000 people", cex=1)   
text(-157.4685, 1.709, "Village with 500 people", cex=1)   
#with(villages[!villages$village=="legend",], text(lon, lat, label=village, pos=1, offset=0.5, font=2, col="black"))
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=3.0)
#with(sites_platy, text(lon, lat, label=site, cex=0.7))
legend(-157.23, 2.08,legend=levels(sites$f.pressure), pt.bg=c("#8c510a","#d8b365","#c7eae5","#5ab4ac","#01665e"), pch=21, bty="n", pt.cex=2.0, cex=1)
text(-157.24, 2.035, "Human disturbance", srt=90, cex=1)
segments(-157.23, 2.004,-157.23, 2.068)
text(-157.3, 1.9, "Bay\nof\nWrecks")
text(-157.52, 1.82, "Vaskess\nBay")
#source("ki_map_files/KI_base_inset.R")

dev.off()


############### same thing as before but making the text bigger ####################


###village data
villages<-read.csv("ki_map_files/KI_villagesDCC_2015update.csv", header = TRUE) # you get an error but it works

#villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.692, 1.722, 1.75))
#villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.5585, -157.5585, -157.5585)
#villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.665, 1.691, 1.715))  # village population on two lines
#villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.548, -157.548, -157.548)   # village population on two lines
villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.665, 1.691, 1.715))  # just village
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.556, -157.556, -157.556)   # just village
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"

setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")
#tiff(file="KI_map_platysites_villages_bigger.tiff",width = 7.6, height = 7.2,units="in",res=300)
#jpeg(file="KI_map_platysites_villages_bigger.jpeg",width = 7.6, height = 7.2,units="in",res=300)
jpeg(file="KI_map_platysites_villages_bigger_3.jpeg",width = 7.6, height = 7.2,units="in",res=300)  # just 'village'
#jpeg(file="KI_map_platysites_villages_bigger_2.jpeg",width = 7.6, height = 7.2,units="in",res=300) # village population on two lines
#pdf(file="KI_map_platysites_villages_bigger.pdf", width = 7.5, height =7)
source("ki_map_files/KI_base_B&W_bigger.R")

# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("red", 0.4))
## legend for village size
#text(-157.487, 1.69, "1500 people", cex=0.69)   # if cex=0.8 then x = -157.4475
#text(-157.487, 1.72, "1000 people", cex=0.69)   # if cex=0.8 then x = -157.4475
#text(-157.49, 1.75, "500 people", cex=0.69)   # if cex=0.8 then x = -157.453
#text(-157.59, 1.72, "Village population", srt=90, cex=0.6)
#segments(-157.577, 1.6925,-157.577, 1.75)
#text(-157.4765, 1.663, "1500 people", cex=0.69)   # village population on two lines
#text(-157.4765, 1.6885, "1000 people", cex=0.69)   # village population on two lines
#text(-157.4795, 1.714, "500 people", cex=0.69)   # village population on two lines
#text(-157.588, 1.69, "Village", srt=90, cex=0.6)# village population on two lines
#segments(-157.566, 1.663,-157.566, 1.714)  # village population on two lines
#segments(-157.577, 1.6925,-157.577, 1.75)
text(-157.4835, 1.663, "1500 people", cex=0.69)   # just village
text(-157.4835, 1.6885, "1000 people", cex=0.69)   # just village
text(-157.4875, 1.714, "500 people", cex=0.69)   # just village
text(-157.588, 1.69, "Village", srt=90, cex=0.6)# just village
segments(-157.575, 1.663,-157.575, 1.714)  # just village
#with(villages[!villages$village=="legend",], text(lon, lat, label=village, pos=1, offset=0.5, font=2, col="black"))
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4) 
#with(sites_platy, text(lon, lat, label=site, cex=0.7))
legend(-157.265, 2.075,legend=levels(sites$f.pressure), pt.bg=c("#8c510a","#d8b365","#c7eae5","#5ab4ac","#01665e"), pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.275, 2.01, "Human disturbance", srt=90, cex=0.6)
segments(-157.263, 1.968,-157.263, 2.057)
text(-157.3, 1.88, "Bay of\nWrecks", cex = 0.5)
text(-157.53, 1.82, "Vaskess\nBay", cex = 0.5)
text(-157.592, 2.078, "b", font=2, cex=1.25)
#source("ki_map_files/KI_base_inset.R")

dev.off()



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#######doing versions with the inset so the village legends need to be moved#############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

###village data
villages<-read.csv("ki_map_files/KI_villagesDCC_2015update.csv", header = TRUE) # you get an error but it works

villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.90, 1.932, 1.958))
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.23, -157.229, -157.229)
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"



############### inset and smaller text #####################################
setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")
#tiff(file="KI_map_platysites_villages_inset.tiff",width = 7, height = 7,units="in",res=300)
#jpeg(file="KI_map_platysites_villages_inset.jpeg",width = 7, height = 7,units="in",res=300)
#pdf(file="KI_map_platysites_villages_inset.pdf")
source("ki_map_files/KI_base_B&W.R")

# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("red", 0.4))
## legend for village size
text(-157.176, 1.90, "Village with\n1500 people", cex=1)   
text(-157.176, 1.93, "Village with\n1000 people", cex=1)   
text(-157.176, 1.96, "Village with\n500 people", cex=1)   
#with(villages[!villages$village=="legend",], text(lon, lat, label=village, pos=1, offset=0.5, font=2, col="black"))
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=3.0)
#with(sites_platy, text(lon, lat, label=site, cex=0.7))
legend(-157.23, 2.08,legend=levels(sites$f.pressure), pt.bg=c("#8c510a","#d8b365","#c7eae5","#5ab4ac","#01665e"), pch=21, bty="n", pt.cex=2.0, cex=1)
text(-157.24, 2.035, "Human disturbance", srt=90, cex=1)
segments(-157.23, 2.004,-157.23, 2.068)
text(-157.3, 1.9, "Bay\nof\nWrecks")
text(-157.52, 1.82, "Vaskess\nBay")
par(mar=c(3,4,3,0))
source("ki_map_files/KI_base_inset.R")

dev.off()





############### same thing as before but making the text bigger ####################

setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")
###village data
villages<-read.csv("ki_map_files/KI_villagesDCC_2015update.csv", header = TRUE) # you get an error but it works

villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 1.852, 1.8791, 1.9045))
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183, -157.263, -157.263, -157.263)
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"

setwd("/Users/KristinaTietjen/Documents/Git_Hub/KI_Platy/figures")
#tiff(file="KI_map_platysites_villages_inset_bigger.tiff",width = 7.6, height = 7.2,units="in",res=300)
#jpeg(file="KI_map_platysites_villages_inset_bigger.jpeg",width = 7.5, height = 7,units="in",res=300)
pdf(file="KI_map_platysites_villages_inset_bigger.pdf", width = 7.5, height =7)
source("ki_map_files/KI_base_B&W_bigger.R")

# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("red", 0.4))
## legend for village size
text(-157.19, 1.85, "1500 people", cex=0.66)   
text(-157.19, 1.877, "1000 people", cex=0.66)   
text(-157.1945, 1.904, "500 people", cex=0.66)  
#with(villages[!villages$village=="legend",], text(lon, lat, label=village, pos=1, offset=0.5, font=2, col="black"))
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4) 
#with(sites_platy, text(lon, lat, label=site, cex=0.7))
legend(-157.265, 2.075,legend=levels(sites$f.pressure), pt.bg=c("#8c510a","#d8b365","#c7eae5","#5ab4ac","#01665e"), pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.275, 2.01, "Human disturbance", srt=90, cex=0.6)
segments(-157.263, 1.968,-157.263, 2.057)
text(-157.315, 1.88, "Bay of\nWrecks", cex = 0.45)
text(-157.524, 1.825, "Vaskess\nBay", cex = 0.45)
par(mar=c(1.3,0.9,0.5,0.5))
source("ki_map_files/KI_base_inset_forbigger.R")

dev.off()


