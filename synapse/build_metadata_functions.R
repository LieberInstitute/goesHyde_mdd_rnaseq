
build_metadata <- function(template_xlsx, data, id_col, id) {
  dict <- read_excel(template_xlsx, sheet = "Dictionary")
  dict_id <- dict$key[match(id_col, dict$col)]
  message(paste("data ID:", id_col, "template id:", dict_id))
  # get extract variable cols from data
  data_hasCol <- !is.na(dict$col) & dict$col != "?"
  data_col <- dict$col[data_hasCol]
  
  # message(paste(colnames(data), collapse = ","))
  # message(paste(data_col, collapse = ","))
  # check colnames match
  col_match <- data_col %in% colnames(data)
  if (!all(col_match)) {
    missing <- data_col[!col_match]
    message("Missing cols:", paste(missing, collapse = ","))
    return(NULL)
  } else {
    message(sum(col_match), " matches in data")
  }
  
  # build data Variable
  dataV <- data[data[[id_col]] %in% id, ]
  dataV <- dataV[, data_col]
  colnames(dataV) <- dict$key[data_hasCol]
  dataV <- replace_values(template_xlsx, dataV)
  # build data Same
  dataS <- t(data.frame(dict$value[!data_hasCol]))
  dataS <- do.call("rbind", replicate(length(id), dataS, simplify = FALSE))
  dataS <- cbind(dataS, id)
  colnames(dataS) <- c(dict$key[!data_hasCol], dict_id)
  # build data All
  dataA <- merge(dataV, dataS, by = dict_id)
  temp <- read_excel(template_xlsx, sheet = "Template")
  
  meta_data <- rbind(temp, dataA)
  return(meta_data)
}

get_fastq_info <- function(fastq) {
  l <- system(paste0("zcat ", fastq, ' | grep "@" | head -n 1'), intern = TRUE) %>%
    strsplit(" ") %>%
    unlist()
  
  l1 <- strsplit(l, ":") %>% unlist()
  info_names <- c(
    "instrument", "rna_id", "flow_cell", "flowcell_lane",
    "title_number", "x_cord", "y_cord", "pair", "filtered",
    "control_bits", "index_seq"
  )
  names(l1) <- info_names
  return(l1)
}

replace_value <- function(value_row, dataV) {
  cn <- value_row[1]
  v <- value_row[2]
  lv <- value_row[3]
  dataV[[cn]][dataV[[cn]] == lv] <- v
  return(dataV)
}

make_value_df <- function(template_xlsx) {
  value_df <- read_excel(template_xlsx, "Values") %>%
    select(key, value, LIBD_value) %>%
    filter(
      !is.na(LIBD_value),
      value != LIBD_value
    ) %>%
    as.data.frame()
  message(nrow(value_df), " values to replace")
  return(value_df)
}

replace_values <- function(template_xlsx, dataV) {
  value_df <- make_value_df(template_xlsx)
  nr <- nrow(value_df)
  if (nr == 0) {
    return(dataV)
  } else {
    for (row in 1:nrow(value_df)) {
      dataV <- replace_value(unlist(value_df[row, ]), dataV)
    }
    return(dataV)
  }
}
