library(ggplot2)
library(dplyr)
library(ggrepel)

pd_mdd <- read.csv("../data/GoesMDD_pd_n1146.csv")
mds <- read.table("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/merged_batches_topmed/LIBD_Brain_merged_topmed_042820.mds", header = TRUE)
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
brain_sentrix <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/brain_sentrix_swap.csv") 

mds2 <- mds %>% mutate(ID = paste(FID,IID,sep = "_")) %>%
  left_join(brain_sentrix) %>%
  left_join(lims %>% select(BrNum, Race)) %>%
  mutate(mdd = ID %in% pd_mdd$genoSample,
         mdd_nw = mdd & Race != "CAUC")%>%
  arrange(-mdd_nw) %>% 
  filter(!(is.na(Race) & mdd))

pca_race <- mds2 %>%
  ggplot(aes(C1, C2, color = Race)) +
  facet_wrap(~mdd)+
  geom_point(shape = 1)

ggsave(pca_race, filename = "pca_race.png")

pca_race_mdd_nw <- mds2 %>%
  mutate(l = ifelse(mdd_nw, BrNum, NA))%>%
  ggplot(aes(C1, C2, color = Race, label = l)) +
  facet_wrap(~mdd_nw)+
  geom_text_repel()+
  geom_point(shape = 1)

ggsave(pca_race_mdd_nw, filename = "pca_race_mdd_nw.png")
