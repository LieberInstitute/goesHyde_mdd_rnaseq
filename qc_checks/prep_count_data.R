###
library(SummarizedExperiment)
library(jaffelab)
library(readxl)
library(RColorBrewer)
library(janitor)

## make directory
dir.create("count_data")
dir.create("pdfs")

## load rses
load("preprocessed_data/rse_gene_human_mouse_cross_species_n40.Rdata")

# read in pheno
pd = read_excel("MAC21 RNAseq sample information 2020_03_17.xlsx")
pd = clean_names(pd)
pd = as.data.frame(pd)
colnames(pd)[1] = "SampleID"

# drop first column
mm2 = match(colnames(rse_gene), pd$SampleID)
pheno = cbind(pd[mm2,], colData(rse_gene))
rownames(pheno) = pheno$SAMPLE_ID

### append
colData(rse_gene) = pheno

## relevel
rse_gene$Group = factor(rse_gene$sample_info)
rse_gene$Tissue = factor(gsub(" Fat", "", rse_gene$sample_source_organ))
rse_gene$Tissue = relevel(rse_gene$Tissue, "Liver")

## metrics
rse_gene$hg38_fraction = rse_gene$h_totalMapped/rse_gene$combined_totalMapped
rse_gene$mm10_fraction = rse_gene$m_totalMapped/rse_gene$combined_totalMapped

## gene assignment
sIndexes =splitit(rowData(rse_gene)$gene_species)
speciesSums = sapply(sIndexes, function(ii)  colSums(assays(rse_gene)$count[ii,]))
rse_gene$h_geneRate = speciesSums[,"human"]/rse_gene$h_totalMapped
rse_gene$m_geneRate = speciesSums[,"mouse"]/rse_gene$m_totalMapped

## make group for plots
rse_gene$Label = paste0(rse_gene$Tissue, "_", rse_gene$Group)
rse_gene$Label = factor(rse_gene$Label, levels = c(
	"Liver_Euploid", "Liver_MAC21", "Brown_Euploid", "Brown_MAC21", 
	"SubQ_Euploid", "SubQ_MAC21", "Visceral_Euploid", "Visceral_MAC21"))

#### QC ####
pdf("pdfs/alignFraction_vs_label.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(13,6,2,2), cex.axis=1.8,cex.lab=1.8)
boxplot(rse_gene$mm10_fraction ~ rse_gene$Label, las=3,xlab="",
	outline = FALSE, ylab = "mouse map fraction")
points(rse_gene$mm10_fraction ~ jitter(as.numeric(rse_gene$Label),amount=0.15),
	pch = 21, bg= rse_gene$Group)
boxplot(rse_gene$hg38_fraction ~ rse_gene$Label, las=3,xlab="",
	outline = FALSE, ylab = "human map fraction")
points(rse_gene$hg38_fraction ~ jitter(as.numeric(rse_gene$Label),amount=0.15),
	pch = 21, bg= rse_gene$Group)
dev.off()

pdf("pdfs/assignFraction_vs_label.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(13,6,2,2), cex.axis=1.8,cex.lab=1.8)
boxplot(rse_gene$m_geneRate ~ rse_gene$Label, las=3,xlab="",
	outline = FALSE, ylab = "mouse assign fraction")
points(rse_gene$m_geneRate ~ jitter(as.numeric(rse_gene$Label),amount=0.15),
	pch = 21, bg= rse_gene$Group)
boxplot(rse_gene$h_geneRate ~ rse_gene$Label, las=3,xlab="",
	outline = FALSE, ylab = "human assign fraction")
points(rse_gene$h_geneRate ~ jitter(as.numeric(rse_gene$Label),amount=0.15),
	pch = 21, bg= rse_gene$Group)
dev.off()

### save
save(rse_gene, file="count_data/wong_MAC21_rse_gene_annotated_n40.Rdata")


#######################
## update exons ###
load("preprocessed_data/rse_exon_human_mouse_cross_species_n40.Rdata")
rse_exon = rse_exon[,colnames(rse_gene)]
colData(rse_exon) = colData(rse_gene)
save(rse_exon, file="count_data/wong_MAC21_rse_exon_annotated_n40.Rdata")
