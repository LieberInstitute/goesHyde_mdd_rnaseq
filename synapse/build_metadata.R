library("readxl")
library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")
library("dplyr")
library("reshape2")
library("purrr")

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
  dataV <- replace_values(template_xlsx,dataV)
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

get_fastq_info <- function(fastq){
  l = system(paste0('zcat ', fastq ,' | grep "@" | head -n 1'), intern=TRUE) %>%
    strsplit(" ") %>%
    unlist()
  
  l1 <- strsplit(l,":") %>% unlist
  info_names <- c("instrument","rna_id","flow_cell","flowcell_lane",
                  "title_number","x_cord","y_cord","pair","filtered",
                  "control_bits","index_seq")
  names(l1) <- info_names
  return(l1)
}

replace_value <- function(value_row, dataV){
  cn <- value_row[1]
  v <- value_row[2]
  lv <- value_row[3]
  dataV[[cn]][dataV[[cn]] == lv] <- v
  return(dataV)
}

make_value_df <- function(template_xlsx){
  value_df <- read_excel(template_xlsx,"Values") %>%
    select(key,value, LIBD_value) %>%
    filter(!is.na(LIBD_value),
           value != LIBD_value) %>%
    as.data.frame()
  message(nrow(value_df)," values to replace")
  return(value_df)}

replace_values <- function(template_xlsx, dataV){
  value_df <- make_value_df(template_xlsx)
  nr <- nrow(value_df)
  if (nr == 0){
    return(dataV)
  }else{
      for (row in 1:nrow(value_df)) {
        dataV <- replace_value(unlist(value_df[row,]),dataV)
      }
      return(dataV)
    }
}

#### load data ####
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv", as.is=TRUE)

# add brain weight to pd
brain_weight <- lims %>% select(BrNum, `Brain.Weight..gram.`)

pd_mdd <- read.csv("../data/GoesMDD_pd_n1140.csv") %>% 
  filter(Experiment == "psychENCODE_MDD") %>%
  left_join(brain_weight, by = "BrNum")  %>%
  mutate(BrodmannArea =  ifelse(BrainRegion == "anterior cingulate cortex",25,NA))

# add imputed Sex to lims
brain_sentrix <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/brain_sentrix_swap.csv") %>%
  filter(ID %in% pd_mdd$genoSample)

lims <- lims %>% left_join(brain_sentrix %>% select(BrNum, genoSex)) 

# replace values 
gia <- read.csv("genotypeInferredAncestry.csv")

lims <- lims %>%
  left_join(gia) %>%
  mutate(genotypeInferredAncestry = ifelse(is.na(genotypeInferredAncestry),Race,genotypeInferredAncestry)) %>%
  select(-PrimaryDx) %>%
  inner_join(pd_mdd %>% select(PrimaryDx, BrNum) %>% unique, by = "BrNum")

#build manifest
mdd_manifest <- pd_mdd %>%
  select(RNum, BrNum, read1, read2) %>%
  melt(id.vars = c("RNum", "BrNum")) %>%
  select(RNum, BrNum, path = value) %>%
  mutate(metadataType = NA, fileFormat = "fastq") %>%
  filter(!is.na(path))

#fastq info
fastq_info <- map(mdd_manifest$path, get_fastq_info)
fastq_info_df <- fastq_info %>% as.data.frame() %>% t()
fastq_info_df <- cbind(fastq_info_df,mdd_manifest[,"RNum", drop=FALSE])
rownames(fastq_info_df) <- NULL

flow_cell <- fastq_info_df %>% select(flow_cell, RNum) %>% unique
pd <- pd  %>% left_join(flow_cell, by = "RNum") %>%
  filter(Dataset == "psychENCODE_MDD")
# list samples
BrNum_mdd <- unique(pd_mdd$BrNum)
RNum_mdd <- pd_mdd$RNum

#### Build metadata ####
# Individual
indi_md <- build_metadata("template_individual_human.xlsx", lims, "BrNum", BrNum_mdd)
dim(indi_md)
# [1] 317  33
indi_md_fn <- "psychENCODE_MDD_individual_human.csv"
write.csv(indi_md, file = indi_md_fn, row.names = FALSE)

# Biospecimen
bio_md <- build_metadata("template_biospecimen.xlsx", pd_mdd, "RNum", RNum_mdd)
dim(bio_md)                                                                                                                                                                                      
# [1] 627  12
bio_md_fn <- "psychENCODE_MDD_biospecimen.csv"
write.csv(bio_md, file = bio_md_fn, row.names = FALSE)

# Assay
assay_md <- build_metadata("template_assay_rnaSeq.xlsx", pd, "RNum", RNum_mdd)
dim(assay_md)
# [1] 627  19
assay_md_fn <- "psychENCODE_MDD_assay_rnaSeq.csv"
write.csv(assay_md, file = assay_md_fn, row.names = FALSE)

# Manifest
mani_md_fn <- "psychENCODE_MDD_manifest.tsv"

meta_files <- c(indi_md_fn, bio_md_fn, assay_md_fn, mani_md_fn)
meta_paths <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/synapse/", meta_files)
meta_manifest <- data.frame(rep(NA,4),
                            rep(NA,4),
                            meta_paths, 
                            c("individual", "biospecimen", "assay", "manifest"),
                            ss(meta_files,"\\.",2))
colnames(meta_manifest) <- c("RNum","BrNum","path", "metadataType","fileFormat")

mdd_all_mani <- rbind(meta_manifest, mdd_manifest) %>%
  mutate(assay = ifelse(fileFormat == "fastq","rnaSeq",NA))
mani_md <- build_metadata("template_manifest.xlsx", mdd_all_mani, "path", mdd_all_mani$path)
dim(mani_md)
# [1] 1258   16
write.table(mani_md, file = mani_md_fn, row.names=FALSE, sep="\t")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('build_metadata', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript build_metadata.R")
