library("readxl")
library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")
library("dplyr")

lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd <- read.csv("../data/GoesMDD_pd_n1146.csv")


build_metadata <- function(template_xlsx, data, id_col, id){
  dict <- read_excel(template_xlsx, sheet = "Dictionary")
  dict_id <- dict$key[match(id_col, dict$col)]
  message(paste("data ID:", id_col, "template id:",dict_id))
  #get extract variable cols from data
  data_hasCol <- !is.na(dict$col) & dict$col != "?"
  data_col <- dict$col[data_hasCol]
  message(all(data_col %in% colnames(data)))
  
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

indi_md <- build_metadata("template_individual_human.xlsx", lims, "BrNum", indi_ID)

