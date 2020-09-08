####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(segmented)

## load
load("../data/rse_gene_GoesZandi_n1100.rda")
load("../data/rse_exon_GoesZandi_n1100.rda")
load("../data/rse_jxn_GoesZandi_n1100.rda")
load("../data/rse_tx_GoesZandi_n1100.rda")


####################
### Regions combined
exprs <- list(
    'Gene' = recount::getRPKM(rse_gene, "Length"),
    'Exon' = recount::getRPKM(rse_exon, "Length"),
    'Jxn' = recount::getRPKM(rse_jxn, "Length"),	# already set to 100
    'Tx' = assays(rse_tx)$tpm
)

## Identify potential cutoffs
seed <- 20191217
seeds <- seed + 0:3
names(seeds) <- names(exprs)
cutoffs <- sapply(names(exprs), function(type) {
    message(type)
    # pdf(paste0('suggested_expr_cutoffs_', tolower(type), '.pdf'), width = 12)
    cuts <- jaffelab::expression_cutoff(exprs[[type]], seed = seeds[type])
    message(paste(cuts, collapse = ' '))
    cut <- max(cuts)
    # dev.off()
    return(cut)
})
# Gene
# 2019-12-17 16:22:30 the suggested expression cutoff is 0.2
# 0.25 0.16
# Exon
# 2019-12-17 16:23:54 the suggested expression cutoff is 0.24
# 0.29 0.2
# Jxn
# 2019-12-17 16:25:37 the suggested expression cutoff is 0.3
# 0.22 0.37
# Tx
# 2019-12-17 16:26:22 the suggested expression cutoff is 0.3
# 0.38 0.22

cutoffs
# Gene Exon  Jxn   Tx
# 0.25 0.29 0.37 0.38



### Filter RSEs
means <- lapply(exprs, rowMeans)

rowRanges(rse_gene)$meanExprs <- means[['Gene']]
rowRanges(rse_gene)$passExprsCut <- means[['Gene']] > cutoffs['Gene']
rse_gene <- rse_gene[rowRanges(rse_gene)$passExprsCut]
save(rse_gene, file = 'rse_gene.Rdata')

rowRanges(rse_exon)$meanExprs <- means[['Exon']]
rowRanges(rse_exon)$passExprsCut <- means[['Exon']] > cutoffs['Exon']
rse_exon <- rse_exon[rowRanges(rse_exon)$passExprsCut]
save(rse_exon, file = 'rse_exon.Rdata')

rowRanges(rse_jxn)$meanExprs <- means[['Jxn']]
rowRanges(rse_jxn)$passExprsCut <- means[['Jxn']] > cutoffs['Jxn']
rse_jxn <- rse_jxn[rowRanges(rse_jxn)$passExprsCut]
save(rse_jxn, file = 'rse_jxn.Rdata')

rowRanges(rse_tx)$meanExprs <- means[['Tx']]
rowRanges(rse_tx)$passExprsCut <- means[['Tx']] > cutoffs['Tx']
rse_tx <- rse_tx[rowRanges(rse_tx)$passExprsCut]
save(rse_tx, file = 'rse_tx.Rdata')



# sgejobs::job_single('get_expression_cutoffs', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript get_expression_cutoffs.R")







