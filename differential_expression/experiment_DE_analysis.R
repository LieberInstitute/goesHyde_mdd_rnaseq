
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
                     exp_t = outGene_ol$t, exp_log10_p = -log10(outGene_ol$P.Value))

t_sacc <- data.frame(de_F = outGene_sACC$F, de_log10_p = -log10(outGene_sACC$P.Value),
                     exp_t = outGene_ol$t, exp_log10_p = -log10(outGene_ol$P.Value))

t_amyg_scatter <- ggplot(t_amyg ,aes(x = de_F, y = exp_t)) +
  geom_point(size = .5) +
  labs(title = "DE vs. Experiment", subtitle = "Amyg",
       x = "Full DE F-stat",
       y = "Experimental DE t-stat")
ggsave(t_amyg_scatter, filename = "plots/DE_vs_Experiment_Amyg.png")

log10p_amyg_scatter <- ggplot(t_amyg ,aes(x = de_log10_p, y = exp_log10_p)) +
  geom_point(size = .5) +
  labs(title = "DE vs. Experiment", subtitle = "Amyg",
       x = "Full DE -log10(p)",
       y = "Experimental DE -log10(p)")
ggsave(log10p_amyg_scatter, filename = "plots/log10p_DE_vs_Experiment_Amyg.png")

t_sacc_scatter <- ggplot(t_sacc ,aes(x = de_F, y = exp_t)) +
  geom_point(size = .5) +
  labs(title = "DE vs. Experiment", subtitle = "sACC",
       x = "Full DE F-stat",
       y = "Experimental DE t-stat")
ggsave(t_sacc_scatter, filename = "plots/DE_vs_Experiment_sACC.png")

log10p_sacc_scatter <- ggplot(t_sacc ,aes(x = de_log10_p, y = exp_log10_p)) +
  geom_point(size = .5) +
  labs(title = "DE vs. Experiment", subtitle = "Amyg",
       x = "Full DE -log10(p)",
       y = "Experimental DE -log10(p)")
ggsave(log10p_sacc_scatter, filename = "plots/log10p_DE_vs_Experiment_sACC.png")

#sgejobs::job_single('experiment_DE_analysis', create_shell = TRUE, memory = '5G', command = "Rscript experiment_DE_analysis.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
