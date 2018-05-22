load("data/KI_Compartment_f_coral_grouped.RData")
sample <- sample_data(phy97.f.c)
write.csv(x=sample,file="data/sample_metadata.csv")

sample_coral <- sample_data(phy97.f.c.coral)
write.csv(x=sample_coral,file="data/sample_coral_metadata.csv")
