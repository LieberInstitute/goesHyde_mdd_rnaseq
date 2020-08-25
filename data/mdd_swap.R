library(jaffelab)
library(VariantAnnotation)
library(readxl)
library(Rsamtools)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(session_info)
#### Load Data ####
load("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/corLong2.Rdata", verbose = TRUE)
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv", as.is=TRUE)

# MDD data
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/rse_gene_raw_GoesZandi_n1174.rda", verbose = TRUE)
mdd_pd <- colData(rse_gene)
dim(mdd_pd) #[1] 1174   13
length(unique(mdd_pd$BrNum)) #[1] 607

table(rse_gene$Experiment)
# GoesMDD ZandiBPD 
# 634      540

#### Swap and drop rna samples ###
# load pd data
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv") %>%
  rename(AgeDeath = Age, BrainRegion = Region, PrimaryDx = Dx, Experiment = Dataset) %>%
  select(colnames(mdd_pd)) %>%
  filter(RNum %in% mdd_pd$RNum) 

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
nrow(brain_sentrix1) ==length(unique(brain_sentrix1$BrNum)) 

#### Add brain sentrix ID to mdd_data ####
pd <- pd %>% left_join(brain_sentrix1, by = "BrNum")

# check for good cor in corLong2
corLong2_mdd <- corLong2 %>% select(RNum, genoSample, cor) %>%
  filter(genoSample %in% pd$genoSample,
         RNum %in% pd$RNum) %>%
  group_by(RNum,genoSample) %>%
  arrange(-cor)%>%
  slice(1) %>%
  ungroup

pd_check <- pd %>% left_join(corLong2_mdd, by = c("RNum","genoSample"))

table(pd_check$cor >= 0.59 & !is.na(pd_check$cor))
#Drop 18 samples 
# FALSE  TRUE 
# 18  1146 

# drop break down by Experiment
pd_check %>%
  count(Experiment, dna_match = cor >= 0.59 & !is.na(cor), good_rna = BrNum != "drop")
# Experiment dna_match good_rna     n
# <chr>      <lgl>     <lgl>    <int>
# 1 GoesMDD    FALSE     FALSE        3
# 2 GoesMDD    FALSE     TRUE         4
# 3 GoesMDD    TRUE      TRUE       627
# 4 ZandiBPD   FALSE     FALSE        7
# 5 ZandiBPD   FALSE     TRUE         4
# 6 ZandiBPD   TRUE      TRUE       519

pd_good <- pd_check %>% 
  filter(BrNum != "drop" & cor >=0.59) 

#filter out brains not in lims
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
br_missing_lims <- pd_good$BrNum[!pd_good$BrNum %in% lims$BrNum] %>% unique
message(paste("Brains missing from lims:",length(br_missing_lims)))
message(paste(br_missing_lims, collapse = " ,"))

pd_good <- pd_good %>% filter(BrNum %in% lims$BrNum)

# create final pd table
pd <- pd_good %>% select(-cor)

dna_rna_cor_histo <- pd_good %>%
  ggplot(aes(cor)) +
  geom_histogram() +
  labs(title = "DNA vs. RNA cor",
       subtitle = paste0(nrow(pd)," good samples for MDD"))

ggsave(filename = "MDDsamples_cor_histo.jpg", plot = dna_rna_cor_histo)



#MDD fastq files
fastq_mdd <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],".fastq.gz")
fastq_mdd2 <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],"_read2.fastq.gz")
#check valid paths
fe_mdd <- file.exists(c(fastq_mdd, fastq_mdd2))
table(fe_mdd)
# TRUE 
# 1254 

fastd_df <- data.frame(pd$RNum[pd$Experiment == "psychENCODE_MDD"], fastq_mdd, fastq_mdd2)
colnames(fastd_df) <- c("RNum", "read1", "read2")

# Add fastq paths 
pd <- pd %>% 
  left_join(fastd_df) %>%
  arrange(BrNum) %>%
  arrange(Experiment)

n = nrow(pd) 

fn = paste0("GoesMDD_pd_n",n,".csv")
write.csv(pd, fn)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('mdd_swap', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript mdd_swap.R")

