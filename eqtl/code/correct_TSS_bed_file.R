library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(sessioninfo)
library(here)
library(tidyverse)

options("width"=200)


load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

row_rse <- as.data.frame(rowRanges(rse_gene))
names(row_rse)
row_rse$TSS <- ifelse(row_rse$strand == "+", row_rse$start, row_rse$end)


row_rse <- row_rse[, c(ncol(row_rse), 1:(ncol(row_rse)-1))]

path <- "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/"

gene_sacc = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/gene_sACC.bed.gz"
gene_sacc <- read.table(gzfile(gene_sacc), header = F)  
colnames(gene_sacc)[4] <- "gencodeID"
colnames(gene_sacc)[1] <- "chr"
corner(gene_sacc)
dim(gene_sacc) #25212 5555

gene_sacc <- left_join(gene_sacc, row_rse[, c("gencodeID", "TSS", "strand")],  by = "gencodeID")
gene_sacc$chr_n <- sub("chr", "", gene_sacc$chr)
gene_sacc <- gene_sacc[, c(556:558, 1:555)]
gene_sacc$V2 = gene_sacc$TSS
gene_sacc$V3 = gene_sacc$V2 +1


gene_sacc <- gene_sacc[order(gene_sacc$chr_n, gene_sacc$TSS), ]
gene_sacc <- gene_sacc[, -c(1:3)]

write.table(gene_sacc, file = paste0(path, "updated_TSS/gene_sACC.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#### linux 

module load htslib

cd /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/updated_TSS/
zcat ../gene_sACC.bed.gz | head -1 > header.txt
cat header.txt gene_sACC.bed | bgzip > gene_sACC.bed.gz
tabix -p bed gene_sACC.bed.gz


#Chr    start   end     ID      4463344452_R01C01       4463344439_R01C02       4463344384_R01C01       4572348357_R01C01       5535522118_R01C02       4463344373_R01C02
chr1    17436   17437   ENSG00000278267.1       6       21      10      2       8       9
chr1    29570   29571   ENSG00000227232.5       67      72      31      50      35      22
chr1    137965  137966  ENSG00000269981.1       2       9       6       0       6       4
chr1    200322  200323  ENSG00000279457.3       99      88      111     54      60      32
chr1    297502  297503  ENSG00000228463.9       131     94      55      72      82      17
chr1    348366  348367  ENSG00000236679.2       15      15      22      6       13      9
chr1    629062  629063  ENSG00000225972.1       5       5       11      7       11      3
chr1    629640  629641  ENSG00000225630.1       699     672     882     806     966     574
chr1    631074  631075  ENSG00000237973.1       736     702     898     1060    1161    474


gencodeGTF = read.table("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")

(svnR-4.2.x) 13:32 eqtl $ cat gencode.v25.annotationGRCh38.gtf | grep 'NEGR1' | grep -w 'gene'
chr1    HAVANA  gene    71395940        72282734        .       -       .       gene_id "ENSG00000172260.14"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "NEGR1"; level 2; tag "overlapping_locus"; havana_gene "OTTHUMG00000009698.5";
chr1    HAVANA  gene    71794232        71837012        .       -       .       gene_id "ENSG00000228853.1"; gene_type "sense_intronic"; gene_status "KNOWN"; gene_name "NEGR1-IT1"; level 2; tag "overlapping_locus"; havana_gene "OTTHUMG00000009732.1";
(svnR-4.2.x) 13:32 eqtl $ cat gencode.v25.annotationGRCh38.gtf | grep 'HTT'  | grep -w 'gene't
(svnR-4.2.x) 13:32 eqtl $ cat gencode.v25.annotationGRCh38.gtf | grep 'HTT'  | grep -w 'gene'
chr4    HAVANA  gene    3063471 3074514 .       -       .       gene_id "ENSG00000251075.1"; gene_type "antisense"; gene_status "KNOWN"; gene_name "HTT-AS"; level 2; havana_gene "OTTHUMG00000159915.4";
chr4    ENSEMBL gene    3063472 3063701 .       +       .       gene_id "ENSG00000278778.1"; gene_type "misc_RNA"; gene_status "KNOWN"; gene_name "HTT-AS1_1"; level 3;
chr4    ENSEMBL gene    3069977 3070151 .       +       .       gene_id "ENSG00000274693.1"; gene_type "misc_RNA"; gene_status "KNOWN"; gene_name "HTT-AS1_2"; level 3;
chr4    HAVANA  gene    3074681 3243959 .       +       .       gene_id "ENSG00000197386.10"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "HTT"; level 2; havana_gene "OTTHUMG00000159916.5";
(svnR-4.2.x) 13:32 eqtl $ 


> row_rse[row_rse$Symbol == "NEGR1", ]
                        TSS seqnames    start      end  width strand Length          gencodeID       ensemblID      gene_type Symbol EntrezID Class NumTx meanExprs passExprsCut
ENSG00000172260.14 72282734     chr1 71395940 72282734 886795      -  13406 ENSG00000172260.14 ENSG00000172260 protein_coding  NEGR1   257194 InGen     5  20.03996         TRUE
> row_rse[row_rse$Symbol == "HTT", ]
                       TSS seqnames   start     end  width strand Length          gencodeID       ensemblID      gene_type Symbol EntrezID Class NumTx meanExprs passExprsCut
ENSG00000197386.10 3074681     chr4 3074681 3243959 169279      +  18128 ENSG00000197386.10 ENSG00000197386 protein_coding    HTT     3064 InGen    13  16.60009         TRUE


library(rtracklayer)
library(GenomicFeatures)

# Specify the path to the GTF file
gtf_file <- "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf"

# Read the GTF file
gtf <- readGFF(gtf_file)
gtf <- import(gtf_file, format = "gtf")

# Extract TSS coordinates
tss <- promoters(gtf, upstream = 1, downstream = 0, base = "start", use.names = TRUE)


# Print TSS coordinates
for (seqname in seqnames(tss)) {
  tss_coords <- start(tss[[seqname]])
  cat(paste(seqname, tss_coords, sep = "\t"), "\n")
  
  negative strand
}


