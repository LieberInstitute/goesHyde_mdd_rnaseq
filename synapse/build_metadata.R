library("readxl")
library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")
library("dplyr")

lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd <- read.csv("../data/GoesMDD_pd_n1146.csv")


# Individual data
indi_dict <- read_excel("template_individual_human.xlsx", sheet = "Dictionary")
indi_ID <- unique(pd$BrNum)
length(indi_ID)
# [1] 605

#get data from Lims
indi_fromCol <- !is.na(indi_dict$col) & indi_dict$col != "?"
indi_limsCol <- indi_dict$col[indi_fromCol]
all(indi_limsCol %in% colnames(lims))
# [1] TRUE
indi_lims <- subset(lims[,indi_limsCol], BrNum %in% indi_ID)
colnames(indi_lims) <- indi_dict$key[indi_fromCol]

#build table with values
indi_value_df <- t(as.data.frame(indi_dict$value[!indi_fromCol]))
indi_value_df <- indi_value_df[rep(seq_len(nrow(indi_value_df)),length(indi_ID)),]                  
colnames(indi_value_df) <- indi_dict$key[!indi_fromCol]
rownames(indi_value_df) <- indi_ID

left_join(indi_lims, indi_value_df, by.x = "individualID", by.y = 0)
