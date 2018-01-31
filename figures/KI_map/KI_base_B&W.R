## Script for creating KI map

library(maptools)
library(maps)
library(scales)
library(RColorBrewer)
# library(mapdata)
# library(mapproj)
# library(ggmap)
# library(png)
# library(grImport)
 library(rgdal)
# library(scales)
# library(PBSmapping)


diva_shp<-readShapePoly("figures/KI_map/shapes/diva-gis/KIR_adm0.shp")


land<-brewer.pal(8, "Greys")[4]


############ BASE MAP FOR ALL PLOTS ##########################
mat<-matrix(1, 100, 100)
mat[60:95, 10:40]<-2
mat

layout(mat)
par(mar=c(5,5.1,0.3,0.3))
plot(diva_shp, xlim=c(-157.59, -157.15), ylim=c(1.73, 2.01), col=land, border=alpha("black",0.5), lwd=1)
par(cex=2.5)

axis(1, padj = 0.8, labels=c(expression(paste(157.6, degree, "W")),expression(paste(157.4, degree, "W")), expression(paste(157.2, degree, "W"))), 
     at=c(-157.6, -157.4,  -157.2))
par(cex=2.5)
axis(2, at=c(1.6, 1.7, 1.8, 1.9, 2, 2.1), las=0, labels=c(expression(paste(1.6, degree, "N")),expression(paste(1.7, degree, "N")),expression(paste(1.8, degree, "N")), expression(paste(1.9, degree, "N")),expression(paste(2, degree, "N")),expression(paste(2.1, degree, "N"))))
box(lwd=1)

# # scale bar. for distances at the equator, where 0.1 degree of longitude = 11.132 km 

### at bottom of plot
rect(-157.29,1.65,  -157.192, 1.65, col="black")  ###  EQUAL TO 10 KM
rect(-157.29, 1.65, -157.29, 1.66)
rect(-157.192, 1.65, -157.192, 1.66)
rect(-157.241, 1.65, -157.241, 1.66)
# rect(-157.3102, 1.69, -157.3102, 1.691)
text(-157.29, 1.67, "0km", cex=0.7, font=1)
text(-157.241, 1.67, "5km", cex=0.7, font=1)
text(-157.192, 1.67, "10km", cex=0.7, font=1)
# text(-157.3102, 1.72, "20 km", cex=0.8, font=1)

# ## at top of plot
# rect(-157.3,2.04,  -157.398, 2.04, col="black")  ###  EQUAL TO 10 KM
# rect(-157.3, 2.04, -157.3, 2.05)
# rect(-157.398, 2.04, -157.398, 2.05)
# # rect(-157.3102, 1.69, -157.3102, 1.691)
# text(-157.398, 2.06, "0km", cex=0.8, font=1)
# text(-157.3, 2.06, "10km", cex=0.8, font=1)
# # text(-157.3102, 1.72, "20 km", cex=0.8, font=1)

# draw polygon for north arrow
polygon(y=c(2.01, 2.05, 2.02), x=c(-157.57, -157.56, -157.56), col="black" )
polygon(y=c(2.05, 2.01, 2.02), x=c(-157.56, -157.55, -157.56), col="black" )
text( -157.56,2.06, "N", cex=0.7)

