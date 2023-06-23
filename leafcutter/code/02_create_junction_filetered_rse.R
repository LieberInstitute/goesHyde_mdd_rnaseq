
library(data.table)
library(SummarizedExperiment)
library(tidyverse)
library(here)
library(purrr)

options("width" = 300)


load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)

row_range <- as.data.frame(rowRanges(rse_jxn))

# to prevent scientific notation in output files eg. chr11	5e+05	500483	.	86	-
options(scipen = 999)

counts <- as.data.frame(assay(rse_jxn, "counts"))
counts$jxn <- row.names(counts)
counts[1:5, 1:5]
counts <- counts %>% separate(jxn, into = c("Chr", "coord", "Strand", "empty"), sep = "[:\\(\\)]") # last row is empty ("after )")
counts <- counts %>% separate(coord, into = c("Start", "End"), sep = "[\\-]") 
counts <- counts %>% mutate(chr_n = sub("^chr", "", Chr))
counts <- counts %>% mutate(Start = as.numeric(Start), End = as.numeric(End)) %>% arrange(chr_n, Start)
counts$point <- "."
counts[1:5, 1090:1097]

c(sum(is.na(counts$Chr)), sum(is.na(counts$Start)), sum(is.na(counts$End)))

new_column_names <- str_split(names(counts), "_") %>% sapply(`[`, 1)
names(counts) <- new_column_names

sample_names <- names(counts)[1:1091]

counts[(counts$Start == 5e+05), c("R13950", "Start", "End")]  ## confirm Start is 500000 not 5 e+05...


#LOOP to create new junction files based on filtered rse_object 

  subset_counts <- counts[, c("Chr", "Start", "End", "point", name, "Strand")]
  names(subset_counts)[5] <- "names"
  file_name <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/leafcutter/data/junc_post_processing/", name, ".txt")
  write.table(subset_counts, file = file_name, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Apply the function to each element in sample_names    
walk(sample_names, extract_and_save)
