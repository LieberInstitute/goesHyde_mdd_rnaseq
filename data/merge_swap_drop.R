#####################################

library(edgeR)
library(ggplot2)
library(jaffelab)
#library(readxl)
#library(Rsamtools)
library(sessioninfo)
library(SummarizedExperiment)
#library(sva)
#library(VariantAnnotation)
library(dplyr)

setwd('/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/')

#load objects
load('../preprocessed_data/rse_gene_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_gene					# 634
rse_mdd$Experiment = "GoesMDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE) #changed to 540
rse_bip = rse_gene					# 540
rse_bip$Experiment = "ZandiBPD"

#### MERGE ####
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
rse_bip$BrainRegion = rse_bip$Brain.Region
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
colData(rse_mdd) = colData(rse_mdd)[,colKeep]
colData(rse_bip) = colData(rse_bip)[,colKeep]

# make rowData consistent
rowData(rse_bip)$Symbol = rowData(rse_mdd)$Symbol 	# fill in blank symbols
rowData(rse_bip)$meanExprs = rowData(rse_bip)$gencodeTx = NULL
rowData(rse_mdd)$meanExprs = rowData(rse_mdd)$gencodeTx = NULL

### combine
rse_both = cbind(rse_mdd, rse_bip) #1174

rse_both$PrimaryDx = droplevels(rse_both$PrimaryDx)
rse_both$PrimaryDx = relevel(rse_both$PrimaryDx, ref="MDD")

tempRpkm = recount::getRPKM(rse_both, "Length")
rowData(rse_both)$meanExprs = rowMeans(tempRpkm)

mdd_pd  = colData(rse_both)
table(mdd_pd$BrainRegion,mdd_pd$PrimaryDx)
# rse_gene n1174
rse_gene = rse_both
message("Raw comined data: n" , ncol(rse_gene))


#### Swap and drop rna samples ####
# load pd data
pd_all <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv")

pd <-  pd_all%>%
  rename(AgeDeath = Age, BrainRegion = Region, PrimaryDx = Dx, Experiment = Dataset) %>%
  select(all_of(c(colKeep,"SAMPLE_ID"))) %>%
  filter(RNum %in% mdd_pd$RNum) 

#rownames(pd) <- ss(pd$SAMPLE_ID,";")

pd$Experiment <- factor(pd$Experiment, levels = c("psychENCODE_MDD","psychENCODE_BP"))
table(pd$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 634             540 

# drop 10 overlapping BPD Control samples
pd <- pd %>%
  group_by(RNum) %>% 
  arrange(Experiment) %>%
  slice(1) %>% 
  ungroup()

pd %>% count(Experiment, BrNum == "drop")
# Experiment      `BrNum == "drop"`     n
# <fct>           <lgl>             <int>
#   1 psychENCODE_MDD FALSE               631
# 2 psychENCODE_MDD TRUE                  3
# 3 psychENCODE_BP  FALSE               523
# 4 psychENCODE_BP  TRUE                  7


#### Match with Brain_Sentrix samples ####
#brain sentrix info
brain_sentrix<- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/brain_sentrix_swap.csv") %>%
  filter(BrNum != "drop", BrNum %in% pd$BrNum)

#load batch priority (lower is better)
array_priority <- read.csv("/dcl01/ajaffe/data/lab/brain_swap/DNA_genotyping_array_priorities.csv")
array_priority <- array_priority[,c(1,4)]

# swap filter prioritize
brain_sentrix1 <- brain_sentrix %>%
  left_join(array_priority, by = "Batch") %>%
  group_by(BrNum) %>%
  mutate(nBr = n(),
         n_batch = length(unique(Batch)),
         MinPriority = min(Priority)) %>%
  filter(Priority == MinPriority) %>%
  arrange(ID)%>%
  slice(1) %>%
  ungroup() %>%
  select(genoSample = ID, BrNum, Batch)

# should only have 1 row per BrNum
message("One row per BrNum: ", nrow(brain_sentrix1) ==length(unique(brain_sentrix1$BrNum))) 

#### Add brain sentrix ID to mdd_data ####
pd <- pd %>% left_join(brain_sentrix1, by = "BrNum")

# check for good cor in corLong2
load("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/corLong2.Rdata", verbose = TRUE)

corLong2_mdd <- corLong2 %>% select(RNum, genoSample, cor) %>%
  filter(genoSample %in% pd$genoSample,
         RNum %in% pd$RNum) %>%
  group_by(RNum,genoSample) %>%
  arrange(-cor)%>%
  slice(1) %>%
  ungroup

# check if in lims
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd_check <- pd %>% 
  left_join(corLong2_mdd, by = c("RNum","genoSample")) %>%
  mutate(lims = BrNum %in% lims$BrNum)

# drop break down by Experiment
pd_check %>%
  count(Experiment,
        dna_match = cor >= 0.59 & !is.na(cor), 
        good_rna = BrNum != "drop",
        in_lims = good_rna & lims)
# Experiment      dna_match good_rna in_lims     n
# <fct>           <lgl>     <lgl>    <lgl>   <int>
# 1 psychENCODE_MDD FALSE     FALSE    FALSE       3
# 2 psychENCODE_MDD FALSE     TRUE     TRUE        4
# 3 psychENCODE_MDD TRUE      TRUE     FALSE       4
# 4 psychENCODE_MDD TRUE      TRUE     TRUE      623
# 5 psychENCODE_BP  FALSE     FALSE    FALSE       7
# 6 psychENCODE_BP  FALSE     TRUE     TRUE        4
# 7 psychENCODE_BP  TRUE      TRUE     FALSE       2
# 8 psychENCODE_BP  TRUE      TRUE     TRUE      517

#### DROP ####
pd_good <- pd_check %>% 
  filter(BrNum != "drop" & cor >=0.59 & lims) 
n = nrow(pd_good)
message("After drop: n", n)

# create histogram of rna-dna cor
dna_rna_cor_histo <- pd_good %>%
  ggplot(aes(cor)) +
  geom_histogram() +
  labs(title = "DNA vs. RNA cor",
       subtitle = paste0(nrow(n)," good samples for MDD"))

ggsave(filename = "MDDsamples_cor_histo.jpg", plot = dna_rna_cor_histo)

#Get paths for MDD fastq files
fastq_mdd <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],".fastq.gz")
fastq_mdd2 <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],"_read2.fastq.gz")
#check valid paths
fe_mdd <- file.exists(c(fastq_mdd, fastq_mdd2))
table(fe_mdd)
# TRUE 
# 1246

fastd_df <- data.frame(pd$RNum[pd$Experiment == "psychENCODE_MDD"], fastq_mdd, fastq_mdd2)
colnames(fastd_df) <- c("RNum", "read1", "read2")

# Add fastq paths 
pd <- pd_good %>% 
  left_join(fastd_df) %>%
  arrange(BrNum) %>%
  arrange(Experiment) %>%
  select(-lims)

fn = paste0("GoesMDD_pd_n",n,".csv")
write.csv(pd, fn)

pd <- pd %>% select(all_of(c(colKeep,"genoSample","SAMPLE_ID"))) 
#### Update rse_gene ####
### combine
rse_mdd <- rse_mdd[,which(rse_mdd$RNum %in% pd$RNum)]
rse_bip <- rse_bip[,which(rse_bip$RNum %in% pd$RNum & !rse_bip$RNum %in% rse_mdd$RNum)]
rse_both = cbind(rse_mdd, rse_bip) #1140
pd <- pd[match(rse_both$RNum, pd$RNum),]
ncol(rse_both) == nrow(pd)
colData(rse_both) <- DataFrame(pd)
colnames(rse_both) <- paste0(pd$RNum,"_", pd$Experiment)

rse_gene <- res_both
save(rse_gene, file = paste0("rse_gene_raw_GoesZandi_n",n,".rda"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('merge_swap_drop', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript merge_swap_drop.R")

