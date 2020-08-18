library("readxl")
library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")
library("dplyr")
library("reshape2")


build_metadata <- function(template_xlsx, data, id_col, id){
  dict <- read_excel(template_xlsx, sheet = "Dictionary")
  dict_id <- dict$key[match(id_col, dict$col)]
  message(paste("data ID:", id_col, "template id:",dict_id))
  #get extract variable cols from data
  data_hasCol <- !is.na(dict$col) & dict$col != "?"
  data_col <- dict$col[data_hasCol]
  
  # message(paste(colnames(data), collapse = ","))
  # message(paste(data_col, collapse = ","))
  #check colnames match
  col_match <- data_col %in% colnames(data)
  if(!all(col_match)){
    missing <- data_col[!col_match]
    message("Missing cols:", paste(missing, collapse = ","))
    return(NULL)
  } else {
    message(sum(col_match), " matches in data")
  }
  
  # build data Variable
  dataV <- data[data[[id_col]] %in% id,]
  dataV <- dataV[,data_col]
  colnames(dataV) <- dict$key[data_hasCol]
  
  # build data Same
  dataS <- t(data.frame(dict$value[!data_hasCol]))
  dataS <- do.call("rbind", replicate(length(id), dataS, simplify = FALSE))
  dataS <- cbind(dataS, id)
  colnames(dataS) <-  c(dict$key[!data_hasCol], dict_id)
  #build data All
  dataA <- merge(dataV, dataS, by = dict_id)
  temp <- read_excel(template_xlsx, sheet = "Template")
  
  meta_data <- rbind(temp, dataA)
  return(meta_data)
}

# load data
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
brain_weight <- lims %>% select(BrNum, `Brain.Weight..gram.`)
# ad brain weight to pd
pd_mdd <- read.csv("../data/GoesMDD_pd_n1146.csv") %>% left_join(brain_weight, by = "BrNum")
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv", as.is=TRUE)

BrNum <- unique(pd_mdd$BrNum)
RNum <- pd_mdd$RNum
# Individual
indi_md <- build_metadata("template_individual_human.xlsx", lims, "BrNum", BrNum)
indi_md_fn <- "GoesMDD_individual_human.csv"
write.csv(indi_md, file = indi_md_fn, row.names = FALSE)

# Biospecimen
bio_md <- build_metadata("template_biospecimen.xlsx", pd_mdd, "RNum", RNum)
bio_md_fn <- "GoesMDD_biospecimen.csv"
write.csv(bio_md, file = bio_md_fn, row.names = FALSE)

# Assay
assay_md <- build_metadata("template_assay_rnaSeq.xlsx", pd, "RNum", RNum)
assay_md_fn <- "GoesMDD_assay_rnaSeq.csv"
write.csv(assay_md, file = assay_md_fn, row.names = FALSE)

# Manifest
mani_md_fn <- "GoesMDD_manifest.csv"

mdd_manifest <- pd_mdd %>%
  select(RNum, BrNum, read1, read2) %>%
  melt(id.vars = c("RNum", "BrNum")) %>%
  select(RNum, BrNum, path = value) %>%
  mutate(metadataType = NA) %>%
  filter(!is.na(path))

meta_files <- c(indi_md_fn, bio_md_fn, assay_md_fn, mani_md_fn)
meta_paths <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/synapse", meta_files)
meta_manifest <- data.frame(rep(NA,4),
                            rep(NA,4),
                            meta_paths, 
                            c("individual", "biospecimen", "assay", "manifest"))
colnames(meta_manifest) <- c("RNum","BrNum","path", "metadataType")

mdd_all_mani <- rbind(meta_manifest, mdd_manifest)
mani_md <- build_metadata("template_manifest.xlsx", mdd_all_mani, "path", mdd_all_mani$path)
write.csv(mani_md, file = mani_md_fn, row.names = FALSE)
