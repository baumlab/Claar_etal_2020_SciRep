library(vegan)

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

coral <- vegan_otu(phy97.f.c.coral)
specpool(coral,smallsample = TRUE)
estimateR(coral)
plot(poolaccum(coral))


sediment <- vegan_otu(phy97.f.c.sediment)
specpool(sediment,smallsample = TRUE)
plot(poolaccum(sediment))

water <- vegan_otu(phy97.f.c.water)
specpool(water)
plot(poolaccum(water))


## install iNEXT package from CRAN
install.packages("iNEXT")
library(iNEXT)

i1 <- iNEXT(coral,q = c(0,1,2),datatype = "abundance")
ggiNEXT(i1,type=1,se=TRUE)+ guides(shape=FALSE,color=FALSE,fill=FALSE)+xlim(0,250)+ylim(0,75)

ggiNEXT(i1,type=2,se=TRUE)+ guides(shape=FALSE,color=FALSE,fill=FALSE)+xlim(0,250)+ylim(0,1)

if(!require(iNEXT)) { devtools::install_github("JohnsonHsieh/iNEXT"); library(iNEXT) }


i2<- prepare_inext(otu_table(phy97.f.c.coral))
i3 <- iNEXT(i2,datatype = "incidence_freq")
ggiNEXT(i3,type=1,se=TRUE)+ guides(shape=FALSE,color=FALSE,fill=FALSE)
+xlim(0,500)+ylim(0,5)
+xlim(0,250)+ylim(0,1)
i3$DataInfo

i4<- prepare_inext(otu_table(phy97.f.c.coral.MAeq))
i5 <- iNEXT(i4)

mergedcoral = merge_samples(phy97.f.c.coral, "site")
SD = merge_samples(sample_data(phy97.f.c.coral), "site")
mergedcoral2 <- otu_table(mergedcoral)
i6 <- prepare_inext(mergedcoral2)
i7 <- iNEXT(i6)

mergedcoral3 = merge_samples(phy97.f.c.coral, "Dist")
mergedcoral4 <- otu_table(mergedcoral3)
i8 <- prepare_inext(mergedcoral3)
i9 <- iNEXT(i8)
plot(i9)

mergedMAeq <- merge_samples(phy97.f.c.coral.MAeq,"site")
mergedMAeq2 <- otu_table(mergedMAeq)
i10 <- prepare_inext(mergedMAeq2)
i11<- iNEXT(i10)
plot(i11)

mergedPlob <- merge_samples(phy97.f.c.coral.Plob,"site")
mergedPlob2 <- otu_table(mergedPlob)
i12 <- prepare_inext(mergedPlob2)
i13<- iNEXT(i12)
plot(i13)

mergedPeyd <- merge_samples(phy97.f.c.coral.Peyd,"site")
mergedPeyd2 <- otu_table(mergedPeyd)
i14 <- prepare_inext(mergedPeyd2)
i15<- iNEXT(i14)
plot(i15)

mergedsed <- merge_samples(phy97.f.c.sediment,"site")
mergedsed2 <- otu_table(mergedsed)
i16 <- prepare_inext(mergedsed2)
i17<- iNEXT(i16)
plot(i17)

mergedsed3 <- merge_samples(phy97.f.c.sediment,"Dist")
mergedsed4 <- otu_table(mergedsed3)
i18 <- prepare_inext(mergedsed3)
i19<- iNEXT(i18)
plot(i19)

mergedwat <- merge_samples(phy97.f.c.water,"site")
mergedwat2 <- otu_table(mergedwat)
i20 <- prepare_inext(mergedwat2)
i21<- iNEXT(i20)
plot(i21)

mergedwat3 <- merge_samples(phy97.f.c.water,"Dist")
mergedwat4 <- otu_table(mergedwat3)
i22 <- prepare_inext(mergedwat4)
i23<- iNEXT(i22)
plot(i23)
axis(side=1, at=c(seq(0,20000,1000)))

incid <- as.matrix(colSums(otu_table(phy97.f.c.coral) != 0))
i24 <- iNEXT(incid)
plot(i24)
  
mergedwat5 <- merge_samples(phy97.f.c.water,"field_season")
mergedwat6 <- otu_table(mergedwat5)
i24 <- prepare_inext(mergedwat6)
i25<- iNEXT(i24,q=c(0,1,2))
plot(i25)
