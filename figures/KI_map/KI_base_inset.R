## Script for creating KI map inset

library(maptools)
library(rgdal)
library(scales)

### Complex inset
world_shp<-readShapePoly("figures/KI_map/shapes/ne_110m_land/ne_110m_land")

pdf("figures/KI_map/inset.pdf")

par(mar=c(2.5,4.8,0.5,0.5))

plot(world_shp, xlim=c(-180,-80),ylim=c(-40,40), 
     col=alpha("black", 0.5), bg="white",cex=0.2, lwd=0.2)

box(lwd=1)

axis(1,lwd=0.2, at=c(-180, -140, -100),cex.axis=2, 
     labels=c(expression(paste(180, degree, "W")), 
              expression(paste(140, degree, "W")), 
              expression(paste(100, degree, "W"))))

axis(2, las=2,lwd=0.2, at=c(-40, 0, 40), cex.axis=2, 
     labels=c(expression(paste(40, degree, "S")), 
              expression(paste(0, degree)), 
              expression(paste(40, degree, "N"))))
#rect(-158, 1, -156, 3, border="red", lwd=2)
points(-157.4, 1.9, col=1, pch=2, cex=3)

dev.off()

jpeg("figures/KI_map/inset.jpg",width = 4,height=4,units = "in",res=300)

par(mar=c(2.5,4.5,0.5,0.5))

plot(world_shp, xlim=c(-180,-80),ylim=c(-40,40), 
     col=alpha("black", 0.5), bg="white", cex=0.2, lwd=0.2)

box(lwd=1)

axis(1,lwd=0.2, at=c(-180, -140, -100),cex.axis=2,
     labels=c(expression(paste(180, degree, "W")), 
              expression(paste(140, degree, "W")), 
              expression(paste(100, degree, "W"))))

axis(2, las=2, lwd=0.2, at=c(-40, 0, 40), cex.axis=2, 
     hadj = 0.85, 
     labels=c(expression(paste(40, degree, "S")), 
              expression(paste(0, degree)), 
              expression(paste(40, degree, "N"))))

points(-157.4, 1.9, col=1, pch=2, cex=3)

dev.off()
