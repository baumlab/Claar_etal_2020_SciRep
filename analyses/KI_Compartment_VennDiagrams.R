rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")


# Plot sediment Venn diagrams
venn.diagram(x=list("Before"=sediment.before.types.subclade,"Post-Storm"=sediment.storm.types.subclade,"After"=sediment.after.types.subclade),file = "analyses/Sediment_Venn.jpg", col=timecols[c(1,3,4)],fill=timecols[c(1,3,4)], margin=0,main = "Symbiodinium Types in Sediment",height=6,width=6,units="in",scaled=T)

sed.b <- 6 #length(sediment.before.types.subclade)
sed.s <- 3  #length(sediment.storm.types.subclade)
sed.a <- 6  #length(sediment.after.types.subclade)
sed.bs <- 2 #sum(sediment.before.types.subclade %in% sediment.storm.types.subclade)
sed.sa <- 1 #sum(sediment.storm.types.subclade %in% sediment.after.types.subclade)
sed.ba <- 7 #sum(sediment.before.types.subclade %in% sediment.after.types.subclade)
sed.bsa <- 13 #sum(sediment.before.types.subclade %in% c(sediment.storm.types.subclade, sediment.after.types.subclade))

jpeg(filename="figures/sediment_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag <- euler(c("Before" = (sed.b), "Post-Storm" = (sed.s), "After" = (sed.a), "Before&Post-Storm" = (sed.bs), "Post-Storm&After" = (sed.sa), "Before&After" = (sed.ba), "Before&Post-Storm&After" = (sed.bsa)))
plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
dev.off()


# Plot water Venn diagrams
venn.diagram(x=list("Before"=water.before.types.subclade,"Post-Storm"=water.storm.types.subclade,"After"=water.after.types.subclade),file = "analyses/water_Venn.jpg", col=timecols[c(1,3,4)],fill=timecols[c(1,3,4)], margin=0,main = "Symbiodinium Types in water",height=6,width=6,units="in",scaled=T)

water.b <- 3 #length(water.before.types.subclade)
water.s <- 5  #length(water.storm.types.subclade)
water.a <- 0  #length(water.after.types.subclade)
water.bs <- 0 #sum(water.before.types.subclade %in% water.storm.types.subclade)
water.sa <- 2 #sum(water.storm.types.subclade %in% water.after.types.subclade)
water.ba <- 5  #sum(water.before.types.subclade %in% water.after.types.subclade)
water.bsa <- 11 #sum(water.before.types.subclade %in% c(water.storm.types.subclade, water.after.types.subclade))


jpeg(filename="figures/water_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag2 <- euler(c("Before" = (water.b), "Post-Storm" = (water.s), "After" = (water.a), "Before&Post-Storm" = (water.bs), "Post-Storm&After" = (water.sa), "Before&After" = (water.ba), "Before&Post-Storm&After" = (water.bsa)))
plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
dev.off()



vd1 <- plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = F,main="Sediment")
vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","Post-Storm","After"),main="Water",key=list(space="right",boxes=list(col=timecols), lwd=c(1,4,2),text=list(c("Before","Post-Storm","After"))
     )) # removed: auto.key=T
vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","Post-Storm","After"),main="Water",auto.key = F)

jpeg(filename="figures/water_sediment_venn.jpg", 
     width = 8, height = 4, units="in",res = 300)
grid.arrange(vd1,vd2,ncol=2)
dev.off()


venn.diagram(x=list("Coral"=coral.types.subclade,"Sediment"=sediment.types.subclade,"Water"=water.types.subclade),file = "analyses/coral_sediment_water_Venn.jpg", col=compcols,fill=compcols, margin=0,main = "Symbiodinium Types",height=6,width=6,units="in",scaled=T)

c <- 9
s <- 16
w <- 2
cs <- 7
cw <- 2
sw <- 5
csw <- 17

c+s+w+cs+cw+sw+csw

jpeg(filename="figures/coral_water_sed_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag3 <- euler(c("Coral" = (c), "Sediment" = (s), "Water" = (w), "Coral&Sediment" = (cs), "Coral&Water" = (cw), "Sediment&Water" = (sw), "Coral&Sediment&Water" = (csw)))
plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=compcols,col=compcols,border=compcols,lwd=c(2,2,2))
dev.off()

write.table(file="analyses/coral_subclades.csv",coral.types.subclade,row.names = F,col.names = F,sep=",",quote=F)
write.table(file="analyses/sediment_subclades.csv",sediment.types.subclade,row.names = F,col.names = F,sep=",",quote=F)
write.table(file="analyses/water_subclades.csv",water.types.subclade,row.names = F,col.names = F,sep=",",quote=F)

coral.sediment.subclades <- coral.types.subclade[which((coral.types.subclade %in% sediment.types.subclade) & !(coral.types.subclade %in% water.types.subclade))]
water.sediment.subclades <- water.types.subclade[which((water.types.subclade %in% sediment.types.subclade) & !(water.types.subclade %in% coral.types.subclade))]
coral.water.subclades <- coral.types.subclade[which((coral.types.subclade %in% water.types.subclade) & !(coral.types.subclade %in% sediment.types.subclade))]
coral.only.subclades <- coral.types.subclade[which(!(coral.types.subclade %in% sediment.types.subclade) & !(coral.types.subclade %in% water.types.subclade))]
water.only.subclades <- water.types.subclade[which(!(water.types.subclade %in% sediment.types.subclade) & !(water.types.subclade %in% coral.types.subclade))]
sediment.only.subclades <- sediment.types.subclade[which(!(sediment.types.subclade %in% coral.types.subclade) & !(sediment.types.subclade %in% water.types.subclade))]
coral.sediment.water.subclades <- coral.types.subclade[which((coral.types.subclade %in% sediment.types.subclade) & (coral.types.subclade %in% water.types.subclade))]

