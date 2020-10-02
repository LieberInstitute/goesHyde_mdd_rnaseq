#####################################

library(edgeR)
library(ggplot2)
library(jaffelab)
library(sessioninfo)
library(SummarizedExperiment)
library(dplyr)
library(here)

#load objects
load(here('preprocessed_data','rse_gene_goesHyde_MDD_n634.Rdata'), verbose=TRUE)
rse_mdd = rse_gene					# 634
rse_mdd$Experiment = "psychENCODE_MDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE) #changed to 540
rse_bip = rse_gene					# 540
rse_bip$Experiment = "psychENCODE_BP"

#### MERGE ####
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
rse_bip$BrainRegion = rse_bip$Brain.Region
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment","bamFile")
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

# rse_gene n1174
rse_gene = rse_both
message("Raw comined data: n" , ncol(rse_gene))


#### Swap and flag rna samples to drop####
# load pd data
pd_all <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv")

pd <-  pd_all%>%
  rename(AgeDeath = Age, BrainRegion = Region, PrimaryDx = Dx, Experiment = Dataset) %>%
  select(all_of(c(colKeep,"SAMPLE_ID","rna_preSwap_BrNum"))) %>%
  filter(RNum %in% mdd_pd$RNum) 

pd$Experiment <- factor(pd$Experiment, levels = c("psychENCODE_MDD","psychENCODE_BP"))
table(pd$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 634             540 


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
  select(genoSample = ID, BrNum, Batch, dna_preSwap_BrNum)

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

## Duplicate sample check, identify samples assigned to BrNum\BrainRegion pairs that already exist
dups <- pd %>%
  group_by(BrNum, BrainRegion) %>%
  mutate(nSamp = n()) %>%
  filter(nSamp>1, BrNum != "drop", !is.na(rna_preSwap_BrNum))

# check if in lims and for overlaping controls
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd_check <- pd %>% 
  left_join(corLong2_mdd, by = c("RNum","genoSample")) %>%
  mutate(lims = BrNum %in% lims$BrNum) %>%
  group_by(RNum) %>%
  mutate(overlap = n() >1) %>%
  ungroup

# drop break down by Experiment
pd_check %>%
  mutate(status = ifelse(cor < 0.59 | is.na(cor),"bad dna match",
                         ifelse(BrNum == "drop", "bad rna",
                                ifelse(!lims, "not in lims",
                                       ifelse(RNum %in% dups$RNum ,"swap cauese duplicate",
                                              ifelse(overlap & Experiment == "psychENCODE_BP", 
                                                     "dup control","KEEP")))))) %>%
  group_by(Experiment) %>%
  count(status)
# Experiment      status                    n
# <fct>           <chr>                 <int>
# 1 psychENCODE_MDD bad dna match             7
# 2 psychENCODE_MDD KEEP                    617
# 3 psychENCODE_MDD not in lims               4
# 4 psychENCODE_MDD swap cauese duplicate     6
# 5 psychENCODE_BP  bad dna match            11
# 6 psychENCODE_BP  dup control              10
# 7 psychENCODE_BP  KEEP                    516
# 8 psychENCODE_BP  not in lims               2
# 9 psychENCODE_BP  swap cauese duplicate     1

#### DROP ####
#Get paths for MDD fastq files
fastq_mdd <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],".fastq.gz")
fastq_mdd2 <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],"_read2.fastq.gz")
#check valid paths
fe_mdd <- file.exists(c(fastq_mdd, fastq_mdd2))
message("All ",length(fe_mdd)," fastq files exist:", all(fe_mdd))

fastd_df <- data.frame(pd$RNum[pd$Experiment == "psychENCODE_MDD"], fastq_mdd, fastq_mdd2)
colnames(fastd_df) <- c("RNum", "read1", "read2")

## Add fastq paths 
pd_check <- pd_check %>% 
  left_join(fastd_df, by = "RNum")

## Drop Samples
pd_good <- pd_check %>% 
  filter(BrNum != "drop" & cor >=0.59 & lims & !RNum %in% dups$RNum &
           !(overlap & Experiment == "psychENCODE_BP")) 
n = nrow(pd_good)
message("After drop: n", n)

# create histogram of rna-dna cor
dna_rna_cor_histo <- pd_good %>%
  ggplot(aes(cor)) +
  geom_histogram() +
  labs(title = "DNA vs. RNA cor",
       subtitle = paste0(nrow(n)," good samples for MDD"))

ggsave(filename = "MDDsamples_cor_histo.jpg", plot = dna_rna_cor_histo)

# After drop: n1133
write.csv(pd_good, "raw_GoesZandi_pd.csv")

pd <- pd_good %>% 
  select(all_of(c(colKeep,"genoSample","SAMPLE_ID","rna_preSwap_BrNum","dna_preSwap_BrNum"))) 

#### Create rse_gene for overlap samples ####
pd_overlap <- pd_check %>%
  filter(overlap) %>% 
  select(all_of(c(colKeep,"genoSample","SAMPLE_ID","rna_preSwap_BrNum","dna_preSwap_BrNum"))) %>%
  mutate(colnames = paste0(RNum,"_",Experiment))

table(pd_overlap$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 10              10 

rse_mdd_ol <- rse_mdd[,which(rse_mdd$RNum %in% pd_overlap$RNum)]
rse_bip_ol <- rse_bip[,which(rse_bip$RNum %in% pd_overlap$RNum)]
rse_both_ol = cbind(rse_mdd_ol, rse_bip_ol) #20
colnames(rse_both_ol) <- paste0(rse_both_ol$RNum,"_", rse_both_ol$Experiment)
pd_overlap <- pd_overlap[match(colnames(rse_both_ol),pd_overlap$colnames),]
pd$colnames <- NULL
message("Overlap tables RNum order matches: ",all(rse_both_ol$RNum == pd_overlap$RNum))
message("Overlap tables Experiment order matches: ",all(rse_both_ol$Experiment == pd_overlap$Experiment))
colData(rse_both_ol) <- DataFrame(pd_overlap)

rse_gene <- rse_both_ol
table(rse_gene$Experiment)
save(rse_gene, file = paste0("rse_gene_control_overlap.rda"))

#### Update rse_gene ####
## combine
rse_mdd <- rse_mdd[,which(rse_mdd$RNum %in% pd$RNum)]
rse_bip <- rse_bip[,which(rse_bip$RNum %in% pd$RNum & !rse_bip$RNum %in% rse_mdd$RNum)]
rse_both = cbind(rse_mdd, rse_bip) #1140
pd <- pd[match(rse_both$RNum, pd$RNum),]
ncol(rse_both) == nrow(pd)
colData(rse_both) <- DataFrame(pd)
colnames(rse_both) <- paste0(pd$RNum,"_", pd$Experiment)

rse_gene <- rse_both
save(rse_gene, file = paste0("rse_gene_raw_GoesZandi.rda"))


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('merge_swap_drop', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript merge_swap_drop.R")

