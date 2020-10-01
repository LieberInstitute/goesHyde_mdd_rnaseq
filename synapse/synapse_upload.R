## Instructions here: https://github.com/Sage-Bionetworks/synapserutils#batch-process
install.packages("synapserutils", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
## use transfer node jhpce-transfer01.jhsph.edu

library(synapserutils)
# synLogin("lahuuki")
mani <- read.table("psychENCODE_MDD_manifest.tsv", header = TRUE)
head(mani)

syncToSynapse("psychENCODE_MDD_manifest.tsv", dryRun = TRUE)

syncToSynapse("psychENCODE_MDD_manifest.tsv")
