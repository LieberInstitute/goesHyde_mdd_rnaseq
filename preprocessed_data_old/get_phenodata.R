##
library(jaffelab)

# samples
manifest = read.table("samples.manifest", header=FALSE, stringsAsFactors=FALSE)

# other info
demog = read.csv("Goes_MASTER_file.csv", stringsAsFactors=FALSE)
		
# reorder
pd = demog[match(manifest$V5, demog$RNum),]
pd = pd[,c(3:4,9:12,8,7,1:2,13:18)]
colnames(pd) = gsub("\\.","",colnames(pd))
rownames(pd) = pd$RNum

write.csv(pd, file="samples_phenodata.csv")
