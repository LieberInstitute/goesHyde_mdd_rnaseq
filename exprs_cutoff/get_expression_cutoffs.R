####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(segmented)
library(here)

## load
load(here("data","rse_gene_GoesZandi.rda"), verbose = TRUE)
load(here("data","rse_exon_GoesZandi.rda"), verbose = TRUE)
load(here("data","rse_jxn_GoesZandi.rda"), verbose = TRUE)
load(here("data","rse_tx_GoesZandi.rda"), verbose = TRUE)


#### Add Mds data to rse objects ####
load(here("genotype_data","goesHyde_bipolarMdd_Genotypes_mds.rda"), verbose = TRUE)

## Check that all BrNum are included
message("# BrNum in pd: ", length(unique(pd$BrNum)),
        ", # BrNum in mds: ", nrow(mds))
message("All mds BrNum in pd:", all(rownames(mds) %in% unique(pd$BrNum)))
message("All pd BrNum in mds:", all(unique(pd$BrNum) %in% rownames(mds)))


messge("BrNum order is identical: ", 
       identical(rse_gene$BrNum, rse_exon$BrNum, rse_jxn$BrNum, rse_tx$BrNum))

## reorders according to rse_gene$BrNum  (matched the order of rse_gene)
mds = mds[rse_gene$BrNum,1:5]

## Merge colData
colData(rse_gene) = cbind(colData(rse_gene), mds)
colData(rse_exon) = cbind(colData(rse_exon), mds)
colData(rse_jxn) = cbind(colData(rse_jxn), mds)
colData(rse_tx) = cbind(colData(rse_tx), mds)

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
# 2020-09-08 13:22:57 the suggested expression cutoff is 0.21
# 0.25 0.16
# Exon
# 2020-09-08 13:24:27 the suggested expression cutoff is 0.24
# 0.29 0.2
# Jxn
# 2020-09-08 13:26:21 the suggested expression cutoff is 0.28
# 0.2 0.36
# Tx
# 2020-09-08 13:27:09 the suggested expression cutoff is 0.3
# 0.38 0.22

cutoffs
# Gene Exon  Jxn   Tx 
# 0.25 0.29 0.36 0.38 



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

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()






