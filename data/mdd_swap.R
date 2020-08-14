library(jaffelab)
library(VariantAnnotation)
library(readxl)
library(Rsamtools)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
#### Load Data ####
load("corLong2.Rdata", verbose = TRUE)
pd <- read.csv("pd_swap.csv", as.is=TRUE)
# MDD data
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/rse_gene_raw_GoesZandi_n1174.rda", verbose = TRUE)
mdd_pd <- colData(rse_gene)
dim(mdd_pd) #[1] 1174   13
length(unique(mdd_pd$BrNum)) #[1] 607

#### Swap and drop rna samples ###
pd <- read.csv("pd_swap.csv") %>% select(RNum, BrNum)

mdd_pd <- mdd_pd %>% 
  as.data.frame() %>%
  select(-BrNum) %>%
  left_join(pd, by = c("RNum"))

mdd_pd %>% count(Experiment, BrNum == "drop")
# Experiment `BrNum == "drop"`     n
# <chr>      <lgl>             <int>
#   1 GoesMDD    FALSE               641
# 2 GoesMDD    TRUE                  3
# 3 ZandiBPD   FALSE               543
# 4 ZandiBPD   TRUE                  7


#### Match with Brain_Sentrix samples ####
#brain sentrix info
brain_sentrix<- read.csv("brain_sentrix_swap.csv") %>%
  filter(BrNum != "drop", BrNum %in% mdd_pd$BrNum)

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
mdd_pd <- mdd_pd %>% left_join(brain_sentrix1, by = "BrNum")


# check for good cor in corLong2
corLong2_mdd <- corLong2 %>% select(RNum, genoSample, cor) %>%
  filter(genoSample %in% mdd_pd$genoSample,
         RNum %in% mdd_pd$RNum) %>%
  group_by(RNum,genoSample) %>%
  arrange(-cor)%>%
  slice(1) %>%
  ungroup()

mdd_pd_check <- mdd_pd %>% left_join(corLong2_mdd, by = c("RNum","genoSample"))

table(mdd_pd_check$cor >= 0.59 & !is.na(mdd_pd_check$cor))
#Drop 18 samples 
# FALSE  TRUE 
# 18  1176 

# drop break down by Experiment
mdd_pd_check %>%
  count(cor < 0.59, BrNum == "drop",Experiment)

#how many overlapping 
mdd_pd_check %>% 
  filter(BrNum != "drop" & cor >=0.59) %>%
  count(RNum) %>%
  filter(n>1)

mdd_pd_good <- mdd_pd_check %>% 
  filter(BrNum != "drop" & cor >=0.59) %>%
  group_by(RNum) %>%
  arrange(Experiment) %>%
  slice(1) %>%
  ungroup()

nrow(mdd_pd_good)
mdd_pd <- mdd_pd_good %>% select(-cor)
mdd_pd %>% count(Experiment)
# Experiment     n
# <chr>      <int>
# 1 GoesMDD      627
# 2 ZandiBPD     519

dna_rna_cor_histo <- mdd_pd_good %>%
  ggplot(aes(cor)) +
  geom_histogram() +
  labs(title = "DNA vs. RNA cor",
       subtitle = paste0(nrow(mdd_pd)," good samples for MDD"))

ggsave(filename = "plots/MDDsamples_cor_histo.jpg", plot = dna_rna_cor_histo)



#good
fastq_mdd <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",mdd_pd$RNum[mdd_pd$Experiment == "GoesMDD"],".fastq.gz")
fastq_mdd2 <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",mdd_pd$RNum[mdd_pd$Experiment == "GoesMDD"],"_read2.fastq.gz")
fe_mdd <- file.exists(c(fastq_mdd, fastq_mdd2))
table(fe_mdd)
# TRUE 
# 1254 

fastd_df <- data.frame(mdd_pd$RNum[mdd_pd$Experiment == "GoesMDD"], fastq_mdd, fastq_mdd2)
colnames(fastd_df) <- c("RNum", "read1", "read2")
mdd_pd <- mdd_pd %>% 
  #left_join(fastd_df) %>%
  arrange(BrNum) %>%
  arrange(Experiment)

write.csv(mdd_pd, "GoesMDD_pd_n1146.csv")

summary(mdd_pd$AgeDeath)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.37   34.44   47.33   46.62   55.92   95.27 
