rm(list=ls())

load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")


sediment.before <- subset_samples(phy97.f.c.sediment, data.frame(sample_data(phy97.f.c.sediment))$field_season == "KI2014",prune=TRUE)
sediment.before <- subset_taxa(sediment.before, taxa_sums(sediment.before) > 0, prune=TRUE)

sediment.storm <- subset_samples(phy97.f.c.sediment, data.frame(sample_data(phy97.f.c.sediment))$field_season == "KI2015a_Post", prune=TRUE)
sediment.storm <- subset_taxa(sediment.storm, taxa_sums(sediment.storm) > 0, prune=TRUE)

sediment.after <- subset_samples(phy97.f.c.sediment, data.frame(sample_data(phy97.f.c.sediment))$field_season == "KI2015b", prune=TRUE)
sediment.after <- subset_taxa(sediment.after, taxa_sums(sediment.after) > 0, prune=TRUE)

sediment.before.types <- unique(data.frame(tax_table(sediment.before))$hit)
sediment.storm.types <- unique(data.frame(tax_table(sediment.storm))$hit)
sediment.after.types <- unique(data.frame(tax_table(sediment.after))$hit)

sediment.before.types.subclade <- sediment.before.types
sediment.before.types.subclade <- gsub("_.*","",sediment.before.types.subclade)
sediment.before.types.subclade <- gsub("\\..*","",sediment.before.types.subclade)
sediment.before.types.subclade <- unique(sediment.before.types.subclade)

sediment.storm.types.subclade <- sediment.storm.types
sediment.storm.types.subclade <- gsub("_.*","",sediment.storm.types.subclade)
sediment.storm.types.subclade <- gsub("\\..*","",sediment.storm.types.subclade)
sediment.storm.types.subclade <- unique(sediment.storm.types.subclade)

sediment.after.types.subclade <- sediment.after.types
sediment.after.types.subclade <- gsub("_.*","",sediment.after.types.subclade)
sediment.after.types.subclade <- gsub("\\..*","",sediment.after.types.subclade)
sediment.after.types.subclade <- unique(sediment.after.types.subclade)


venn.diagram(x=list("Before"=sediment.before.types.subclade,"Storm"=sediment.storm.types.subclade,"After"=sediment.after.types.subclade),file = "analyses/Sediment_Venn.jpg", col=timecols[c(1,3,4)],fill=timecols[c(1,3,4)], margin=0,main = "Symbiodinium Types in Sediment",height=6,width=6,units="in",scaled=T)

sed.b <- 6 #length(sediment.before.types.subclade)
sed.s <- 3  #length(sediment.storm.types.subclade)
sed.a <- 6  #length(sediment.after.types.subclade)
sed.bs <- 2 #sum(sediment.before.types.subclade %in% sediment.storm.types.subclade)
sed.sa <- 1 #sum(sediment.storm.types.subclade %in% sediment.after.types.subclade)
sed.ba <- 7 #sum(sediment.before.types.subclade %in% sediment.after.types.subclade)
sed.bsa <- 13 #sum(sediment.before.types.subclade %in% c(sediment.storm.types.subclade, sediment.after.types.subclade))

jpeg(filename="figures/sediment_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag <- euler(c("Before" = (sed.b), "Storm" = (sed.s), "After" = (sed.a), "Before&Storm" = (sed.bs), "Storm&After" = (sed.sa), "Before&After" = (sed.ba), "Before&Storm&After" = (sed.bsa)))
plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
dev.off()

water.before <- subset_samples(phy97.f.c.water, data.frame(sample_data(phy97.f.c.water))$field_season == "KI2014",prune=TRUE)
water.before <- subset_taxa(water.before, taxa_sums(water.before) > 0, prune=TRUE)

water.storm <- subset_samples(phy97.f.c.water, data.frame(sample_data(phy97.f.c.water))$field_season == "KI2015a_Post", prune=TRUE)
water.storm <- subset_taxa(water.storm, taxa_sums(water.storm) > 0, prune=TRUE)

water.after <- subset_samples(phy97.f.c.water, data.frame(sample_data(phy97.f.c.water))$field_season == "KI2015b", prune=TRUE)
water.after <- subset_taxa(water.after, taxa_sums(water.after) > 0, prune=TRUE)

water.before.types <- unique(data.frame(tax_table(water.before))$hit)
water.storm.types <- unique(data.frame(tax_table(water.storm))$hit)
water.after.types <- unique(data.frame(tax_table(water.after))$hit)

water.before.types.subclade <- water.before.types
water.before.types.subclade <- gsub("_.*","",water.before.types.subclade)
water.before.types.subclade <- gsub("\\..*","",water.before.types.subclade)
water.before.types.subclade <- unique(water.before.types.subclade)

water.storm.types.subclade <- water.storm.types
water.storm.types.subclade <- gsub("_.*","",water.storm.types.subclade)
water.storm.types.subclade <- gsub("\\..*","",water.storm.types.subclade)
water.storm.types.subclade <- unique(water.storm.types.subclade)

water.after.types.subclade <- water.after.types
water.after.types.subclade <- gsub("_.*","",water.after.types.subclade)
water.after.types.subclade <- gsub("\\..*","",water.after.types.subclade)
water.after.types.subclade <- unique(water.after.types.subclade)


venn.diagram(x=list("Before"=water.before.types.subclade,"Storm"=water.storm.types.subclade,"After"=water.after.types.subclade),file = "analyses/water_Venn.jpg", col=timecols[c(1,3,4)],fill=timecols[c(1,3,4)], margin=0,main = "Symbiodinium Types in water",height=6,width=6,units="in",scaled=T)

water.b <- 3 #length(water.before.types.subclade)
water.s <- 5  #length(water.storm.types.subclade)
water.a <- 0  #length(water.after.types.subclade)
water.bs <- 0 #sum(water.before.types.subclade %in% water.storm.types.subclade)
water.sa <- 2 #sum(water.storm.types.subclade %in% water.after.types.subclade)
water.ba <- 5  #sum(water.before.types.subclade %in% water.after.types.subclade)
water.bsa <- 11 #sum(water.before.types.subclade %in% c(water.storm.types.subclade, water.after.types.subclade))


jpeg(filename="figures/water_venn.jpg", 
     width = 4, height = 4, units="in",res = 300)
VennDiag2 <- euler(c("Before" = (water.b), "Storm" = (water.s), "After" = (water.a), "Before&Storm" = (water.bs), "Storm&After" = (water.sa), "Before&After" = (water.ba), "Before&Storm&After" = (water.bsa)))
plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = T)
dev.off()



vd1 <- plot(VennDiag, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),auto.key = F,main="Sediment")
vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","During","After"),main="Water",key=list(space="right",boxes=list(col=timecols), lwd=c(1,4,2),text=list(c("Before","During","After"))
     )) # removed: auto.key=T
vd2 <- plot(VennDiag2, quantities = TRUE, font=1, cex=1, alpha=0.5,
            fill=timecols,col=timecols,border=timecols,lwd=c(1,4,2),labels = c("Before","During","After"),main="Water",auto.key = F)

jpeg(filename="figures/water_sediment_venn.jpg", 
     width = 8, height = 4, units="in",res = 300)
grid.arrange(vd1,vd2,ncol=2)
dev.off()


all.types <- unique(data.frame(tax_table(phy97.f.c))$hit)
sediment.types <- unique(data.frame(tax_table(phy97.f.c.sediment))$hit)
water.types <- unique(data.frame(tax_table(phy97.f.c.water))$hit)
coral.types <- unique(data.frame(tax_table(phy97.f.c.coral))$hit)

sediment.types.subclade <- sediment.types
sediment.types.subclade <- gsub("_.*","",sediment.types.subclade)
sediment.types.subclade <- gsub("\\..*","",sediment.types.subclade)
sediment.types.subclade <- unique(sediment.types.subclade)

water.types.subclade <- water.types
water.types.subclade <- gsub("_.*","",water.types.subclade)
water.types.subclade <- gsub("\\..*","",water.types.subclade)
water.types.subclade <- unique(water.types.subclade)

coral.types.subclade <- coral.types
coral.types.subclade <- gsub("_.*","",coral.types.subclade)
coral.types.subclade <- gsub("\\..*","",coral.types.subclade)
coral.types.subclade <- unique(coral.types.subclade)


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

