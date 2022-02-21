library(SummarizedExperiment)
library(tidyverse)
library(here)

load(here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"), verbose = TRUE)
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

#### Check stats for gene 5HTR2a in MDD vs. Control ###

rd <- as.data.frame(rowData(rse_gene))
head(rd)

rd %>%
  filter(grepl("HTR2",Symbol, ignore.case = TRUE))

## Summary Table
pd <- as.data.frame(colData(rse_gene))

table(pd$PrimaryDx)
table(pd$Race)
summary(pd$AgeDeath)

## Summary Table for MDD v Control
mdd_info <- pd %>% 
  filter(PrimaryDx != "Bipolar") %>%
  group_by(BrainRegion, PrimaryDx)  %>%
  summarise(n = n(),
            n_Male = sum(Sex == "M"),
            n_Female = sum(Sex == "F"),
            min_age = min(AgeDeath),
            mean_age = mean(AgeDeath),
            median_age = median(AgeDeath),
            max_age = max(AgeDeath))

write.csv(mdd_info, 
          file = here("differential_expression","data","check_gene","sample_info_MDDvControl.csv"),
          row.names = FALSE)

## Just sep model
outGene <- outGene$sep
## Just MDD
outGene <- transpose(outGene)
outGene <- outGene$MDD
names(outGene)

outGene_combined <- do.call("rbind", map2(outGene, c("Amygdala","sACC"),
                                          ~.x %>% mutate(BrainRegion = .y)))

gene_data <- outGene_combined %>%
  filter(Symbol == "HTR2A") %>%
  as_tibble() %>%
  select(Symbol, gencodeID,BrainRegion, logFC, t, P.Value, FDR = adj.P.Val)
  
write.csv(gene_data, 
          file = here("differential_expression","data","check_gene","HTR2A_MDDvControl.csv"),
          row.names = FALSE)

