

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(sessioninfo)
library(here)
library(dplyr)
library(stringr)
library(rcompanion)

##### Load rse data, examine ####
#to see more column
options("width"=200)

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)


BD <- read.delim("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/prs/mdd_bdp_PGC_BD_PRScs_score.profile", sep = " ")


BD$genoSample <- paste(BD$FID, BD$IID, sep = "_")

#rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]

pd <- as_tibble(colData(rse_gene)) 
pd <- pd %>% distinct(genoSample, .keep_all = TRUE)
str(pd)

dim(BD)

table(BD$genoSample %in% (colData(rse_gene)$genoSample))
table(BD$genoSample %in% (pd$genoSample))

#FALSE  TRUE 
#   21   595 

BD <- BD[(BD$genoSample %in% colData(rse_gene)$genoSample), ]
MDD <- MDD[(MDD$genoSample %in% colData(rse_gene)$genoSample), ]

BD <- BD %>% left_join(pd)
str(BD)
glimpse(BD)
table(BD$PrimaryDx) #total sums 595

BD <- BD[BD$PrimaryDx %in% c("Bipolar", "Control"), ]
BD$PrimaryDx <- droplevels(BD$PrimaryDx)
BD$Dx <- ifelse(BD$PrimaryDx == "Bipolar",1,0)                       
table(BD$Dx)                       
BD <- BD %>% mutate_at(("SCORE"), scale)
base <- glm(Dx ~ snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, family = binomial (link = "logit"), data = BD)
model <- glm(Dx ~ snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + SCORE, family = binomial (link = "logit"), data = BD)
summary(model)  
nagelkerke(model)                      
nagelkerke(base)



##

MDD <- read.delim("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/prs/mdd_bdp_PGC_MDD_PRScs_score.profile", sep = " ")
#to test BD into MDD
#MDD <- read.delim("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/prs/mdd_bdp_PGC_BD_PRScs_score.profile", sep = " ")


MDD$genoSample <- paste(MDD$FID, MDD$IID, sep = "_")
pd <- as_tibble(colData(rse_gene)) 
pd <- pd %>% distinct(genoSample, .keep_all = TRUE)
str(pd)
MDD <- MDD[(MDD$genoSample %in% colData(rse_gene)$genoSample), ]

MDD <- MDD %>% left_join(pd)
str(MDD)
glimpse(MDD)
table(MDD$PrimaryDx) #total sums 595

MDD <- MDD[MDD$PrimaryDx %in% c("MDD", "Control"), ]
MDD$PrimaryDx <- droplevels(MDD$PrimaryDx)
MDD$Dx <- ifelse(MDD$PrimaryDx == "MDD",1,0)                       
table(MDD$Dx)                       
MDD <- MDD %>% mutate_at(("SCORE"), scale)
base <- glm(Dx ~ snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, family = binomial (link = "logit"), data = MDD)
model <- glm(Dx ~ snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + SCORE, family = binomial (link = "logit"), data = MDD)
summary(model)  
nagelkerke(model)                      
nagelkerke(base)

 #0.0671513 - 0.0345596
[1] 0.0325917
> 
    
    
