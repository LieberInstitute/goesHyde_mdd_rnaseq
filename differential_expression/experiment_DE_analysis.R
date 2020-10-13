
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(edgeR)
library(limma)
library(ggplot2)
library(sessioninfo)
library(here)

## Load Data
load("qSVA_MDD_gene_DEresults.rda", verbose = TRUE)

load(here("data","rse_gene_overlap_GoesZandi.rda"), verbose = TRUE)

#### Build Experiment Model ####
mod <- model.matrix(~Experiment + BrNum, data = colData(rse_gene_ol))
colnames(mod)

dge = DGEList(counts = assays(rse_gene_ol)$counts, genes = rowData(rse_gene_ol))
dge <- calcNormFactors(dge)
vGene = voom(dge, mod)
### limma commads
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#compute T-stat
outGene_ol = topTable(eBGene,coef=2,p.value = 1,number=nrow(rse_gene_ol))

table(outGene_ol$adj.P.Val < 0.05)
table(outGene_ol$adj.P.Val < 0.01)

## Save
save(outGene_ol, file = "outGene_overlap.rda")

#### Scatter Plot  ####
## Experiment DE t-stat vs. full Diff Express F-stat
outGene_ol <- outGene_ol[match(outGene_sACC$gencodeID, outGene_ol$gencodeID),]
t_amyg <- data.frame(de_F = outGene_Amyg$F, de_log10_p = -log10(outGene_Amyg$P.Value), 
                     exp_t = outGene_ol$t, exp_log10_p = -log10(outGene_ol$P.Value),
                     BrainRegion = "Amygdyla")

t_sacc <- data.frame(de_F = outGene_sACC$F, de_log10_p = -log10(outGene_sACC$P.Value),
                     exp_t = outGene_ol$t, exp_log10_p = -log10(outGene_ol$P.Value),
                     BrainRegion = "sACC")
t_both <- rbind(t_amyg, t_sacc)

t_sacc_scatter <- ggplot(t_both ,aes(x = de_F, y = exp_t)) +
  geom_point(size = .5) +
  labs(title = "DE Full model vs. Experiment model",
       x = "Full DE F-stat",
       y = "Experimental DE t-stat")+
  facet_wrap(~BrainRegion)
ggsave(t_sacc_scatter, filename = "plots/DE_tstat_vs_Experiment_fstat.png", 
       width = 8, height = 5)

log10p_scatter <- ggplot(t_both ,aes(x = de_log10_p, y = exp_log10_p)) +
  geom_point(size = .5) +
  labs(title = "DE Full model vs. Experiment model",
       x = "Full Data DE -log10(p-value)",
       y = "Experimental DE -log10(p-value)") +
  facet_wrap(~BrainRegion)
ggsave(log10p_scatter, filename = "plots/DE_log10p_vs_Experiment_log10p.png", 
       width = 8, height = 5)

#sgejobs::job_single('experiment_DE_analysis', create_shell = TRUE, memory = '5G', command = "Rscript experiment_DE_analysis.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
