# Clear working environment
rm(list=ls())

# Load necessary packages
library(VennDiagram)
library(eulerr)
library(gridExtra)

# Load necessary data
load("data/KI_Compartment_f_coral_grouped.RData")
load("analyses/KI_Compartment_colors.RData")

# Determine the intersects (i.e. which subclades are shared) between different coral species for Venn plotting
PeydPlob <- intersect(Peyd.ASVs,Plob.ASVs)
PeydPlobMaeq <- intersect(PeydPlob,MAeq.ASVs)
PeydPlob.only <- length(PeydPlob)-length(PeydPlobMaeq)
PeydPlobMaeq.length <- length(PeydPlobMaeq)
PeydMaeq <- intersect(Peyd.ASVs,MAeq.ASVs)
PeydMaeq.only <- length(PeydMaeq)-length(PeydPlobMaeq)
PlobMaeq <- intersect(Plob.ASVs,MAeq.ASVs)
PlobMaeq.only <- length(PlobMaeq)-length(PeydPlobMaeq)
Peyd <- length(Peyd.ASVs) - PeydPlob.only - PeydMaeq.only - PeydPlobMaeq.length
Plob <- length(Plob.ASVs) - PlobMaeq.only - PeydPlob.only - PeydPlobMaeq.length
Maeq <- length(MAeq.ASVs) - PlobMaeq.only - PeydMaeq.only - PeydPlobMaeq.length

## Make figure - jpg
jpeg(filename="figures/Figure_2/Figure_2a_coral_species_venn_OTU.jpg", 
     width = 4, height = 4, units="in",res = 300) # Open jpg
# Create Venn object
VennDiag4 <- euler(c("P. grandis" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. grandis&P. lobata" = (PeydPlob.only), "P. grandis&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
# Plot Venn object
plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off() # Close jpg

## Make figure - pdf
pdf(file="figures/Figure_2/Figure_2a_coral_species_venn_OTU.pdf", 
     width = 4, height = 4) # Open pdf
# Create Venn object
VennDiag4 <- euler(c("P. grandis" = (Peyd), "P. lobata" = (Plob), "M. aequituberculata" = (Maeq), "P. grandis&P. lobata" = (PeydPlob.only), "P. grandis&M. aequituberculata" = (PeydMaeq.only), "P. lobata&M. aequituberculata" = (PlobMaeq.only), "P. grandis&P. lobata&M. aequituberculata" = (PeydPlobMaeq.length)))
# Plot Venn object
plot(VennDiag4, quantities = TRUE, font=1, cex=1, alpha=0.5,
     fill=speccols,col=speccols,border=speccols,lwd=c(2,2,2))
dev.off() # Close pdf
