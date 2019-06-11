rm(list=ls())

library(VennDiagram)
library(eulerr)
library(gridExtra)

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

sed.bs <- intersect(sediment.storm.types.subclade,sediment.before.types.subclade)
sed.bsa <- intersect(sed.bs, sediment.after.types.subclade)
sed.bs.only <- length(sed.bs)-length(sed.bsa)
sed.bsa.length <- length(sed.bsa)
sed.ba <- intersect(sediment.before.types.subclade,sediment.after.types.subclade)
sed.ba.only <- length(sed.ba)-length(sed.bsa)
sed.sa <- intersect(sediment.storm.types.subclade,sediment.after.types.subclade)
sed.sa.only <- length(sed.sa)-length(sed.bsa)
sed.a <- length(sediment.after.types.subclade) - sed.sa.only - sed.ba.only - sed.bsa.length
sed.b <- length(sediment.before.types.subclade) - sed.bs.only - sed.ba.only - sed.bsa.length
sed.s <- length(sediment.storm.types.subclade) - sed.sa.only - sed.bs.only - sed.bsa.length


# jpeg(filename="figures/sediment_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
VennDiag <- euler(c("Before" = (sed.b), "Post-Storm" = (sed.s), "After" = (sed.a), "Before&Post-Storm" = (sed.bs.only), "Post-Storm&After" = (sed.sa.only), "Before&After" = (sed.ba.only), "Before&Post-Storm&After" = (sed.bsa.length)))
plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()


water.bs <- intersect(water.storm.types.subclade,water.before.types.subclade)
water.bsa <- intersect(water.bs, water.after.types.subclade)
water.bs.only <- length(water.bs)-length(water.bsa)
water.bsa.length <- length(water.bsa)
water.ba <- intersect(water.before.types.subclade,water.after.types.subclade)
water.ba.only <- length(water.ba)-length(water.bsa)
water.sa <- intersect(water.storm.types.subclade,water.after.types.subclade)
water.sa.only <- length(water.sa)-length(water.bsa)
water.a <- length(water.after.types.subclade) - water.sa.only - water.ba.only - water.bsa.length
water.b <- length(water.before.types.subclade) - water.bs.only - water.ba.only - water.bsa.length
water.s <- length(water.storm.types.subclade) - water.sa.only - water.bs.only - water.bsa.length



# jpeg(filename="figures/water_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
VennDiag2 <- euler(c("Before" = (water.b), "Post-Storm" = (water.s), "After" = (water.a), "Before&Post-Storm" = (water.bs.only), "Post-Storm&After" = (water.sa.only), "Before&After" = (water.ba.only), "Before&Post-Storm&After" = (water.bsa.length)))
plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()

coral.bs <- intersect(coral.storm.types.subclade,coral.before.types.subclade)
coral.bsa <- intersect(coral.bs, coral.after.types.subclade)
coral.bs.only <- length(coral.bs)-length(coral.bsa)
coral.bsa.length <- length(coral.bsa)
coral.ba <- intersect(coral.before.types.subclade,coral.after.types.subclade)
coral.ba.only <- length(coral.ba)-length(coral.bsa)
coral.sa <- intersect(coral.storm.types.subclade,coral.after.types.subclade)
coral.sa.only <- length(coral.sa)-length(coral.bsa)
coral.a <- length(coral.after.types.subclade) - coral.sa.only - coral.ba.only - coral.bsa.length
coral.b <- length(coral.before.types.subclade) - coral.bs.only - coral.ba.only - coral.bsa.length
coral.s <- length(coral.storm.types.subclade) - coral.sa.only - coral.bs.only - coral.bsa.length


# jpeg(filename="figures/coral_venn.jpg", 
#      width = 4, height = 4, units="in",res = 300)
VennDiag3 <- euler(c("Before" = (coral.b), "Post-Storm" = (coral.s), "After" = (coral.a), "Before&Post-Storm" = (coral.bs.only), "Post-Storm&After" = (coral.sa.only), "Before&After" = (coral.ba.only), "Before&Post-Storm&After" = (coral.bsa.length)))
plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
# dev.off()

vd1 <- plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = F,main="Sediment")
vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","Post-Storm","After"),main="Water",auto.key = F)
vd3 <- plot(VennDiag3, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = F,main="Coral")

class(print(vd1))


jpeg(filename="figures/Figure_5b_sediment_time_venn.jpg",
     width = 4, height = 4, units="in",res = 300)
vd1
dev.off()

jpeg(filename="figures/Figure_5c_water_time_venn.jpg",
     width = 4, height = 4, units="in",res = 300)
vd2
dev.off()

jpeg(filename="figures/Figure_5a_coral_time_venn.jpg",
     width = 4, height = 4, units="in",res = 300)
vd3
dev.off()


pdf(file="figures/Figure_5b_sediment_time_venn.pdf",
     width = 4, height = 4)
vd1
dev.off()

pdf(file="figures/Figure_5c_water_time_venn.pdf",
     width = 4, height = 4)
vd2
dev.off()

pdf(file="figures/Figure_5a_coral_time_venn.pdf",
     width = 4, height = 4)
vd3
dev.off()
# 
# jpeg(filename="figures/Figure_4_venn_by_fieldseason.jpg", 
#      width = 12, height = 4, units="in",res = 300)
# grid.arrange(vd1j)
# dev.off()
# 
