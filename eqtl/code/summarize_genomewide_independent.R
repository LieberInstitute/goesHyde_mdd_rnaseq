library(tidyverse)
library(here)
library(sessioninfo)


#### status update ####
log_files <- list.files(
    path = here("eqtl", "code", "logs"),
    pattern = "tensorqtl_genomewide_independent_argv*",
    full.names = TRUE
)

logs <- map(log_files, readLines)

pair_line <- map_chr(map_chr(logs, 26), ~ gsub("Processing pair: ", "", .x))
names(logs) <- pair_line

last_line <- map_chr(logs, ~ tail(.x, 1))
done <- !grepl("    processing phenotype ", last_line)
message(length(log_files), "/46 logged, ", round(100 * length(log_files) / 46, 1), "%")
message(sum(done), "/46 Done, ", round(100 * sum(done) / 46, 1), "%")

map_chr(last_line[!done], ~ gsub("    processing phenotype ", "", .x))

## 1/24
# 26/43 logged, 60.5%
# 22/43 Done, 51.2%

# Amygdala_chr8 Amygdala_chr10 Amygdala_chr16     sACC_chr16
#     "250/933"       "77/987"     "400/1122"     "730/1110"

## 1/25
# 29/43 logged, 67.4%
# 25/43 Done, 58.1%

# Amygdala_chr8 Amygdala_chr10      sACC_chr7      sACC_chr5
#     "381/933"      "232/987"     "143/1191"      "48/1246"

####
files_inde <- list.files(
    path = here("eqtl", "data", "tensorQTL_out", "genomewide_independent"),
    full.names = TRUE
)

eqtl_cis <- read.csv(here("eqtl", "data", "tensorQTL_out", "genomewide_cis", "gene_Amygdala.csv"))
eqtl_inde <- read.csv(files_inde[[1]], row.names = 1)

head(eqtl_inde)
summary(eqtl_inde$tss_distance)
summary(eqtl_inde$pval_perm)
table(eqtl_inde$pval_perm < 0.01)
# FALSE  TRUE
# 3496    66

#
dim(eqtl_inde)
# [1] 3562   17

## cis =  one result per gene
length(unique(eqtl_inde$phenotype_id))
# [1] 521

eqtl_inde %>%
    count(phenotype_id) %>%
    count(n)

length(unique(eqtl_inde$variant_id))
# [1] 2456

eqtl_inde %>%
    count(variant_id) %>%
    count(n)

# map_int(eqtl_inde, ~length(unique(.x$gencodeID)))
# amyg  sacc
# 25085 25085

# map_int(eqtl_tensor, ~length(unique(.x$variant_id)))
# amyg  sacc
# 16339 16027

test_inde <- head(eqtl_inde)

eqtl_cis %>% filter(variant_id %in% test_inde$variant_id)
