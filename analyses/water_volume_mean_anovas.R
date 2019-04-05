sv <- read.csv("~/../Desktop/sample_vol.csv")

head(sv)

vol_mean <- mean(sv$Volume..L.)
vol_sd <- sd(sv$Volume..L.)
vol_se <- sd(sv$Volume..L.)/sqrt(nrow(sv))

# Compute the analysis of variance
res.aov <- aov(Volume..L. ~ Site.Name, data = sv)
# Summary of the analysis
summary(res.aov)

library(ggplot2)

ggplot(sv) + geom_boxplot(aes(x=Site.Name, y=Volume..L.))

load("~/Data_Analysis/KI_Compartment/data/KI_Compartment_f_coral_grouped.RData")
phy97.f.c.water

ssums <- as.data.frame(sample_sums(phy97.f.c.water))
sitenum <- as.data.frame(sample_data(phy97.f.c.water)$site)
rn <- rownames(sample_data(phy97.f.c.water))
rownames(sitenum) <- rn
ssums_site <- cbind(ssums,sitenum)
colnames(ssums_site) <- c("ssum","site")


# Compute the analysis of variance
res.aov <- aov(ssum ~ site, data = ssums_site)
# Summary of the analysis
summary(res.aov)
