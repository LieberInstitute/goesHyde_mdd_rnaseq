library("tidyverse")
library("jaffelab")
library("sessioninfo")
library("here")

#### Check nums from log files ####
nominal_log_fn <- list.files(here("leafcutter","code", "logs"), pattern = "07_run_tensorqtl", full.names = TRUE)

nominal_log <- map(nominal_log_fn, readLines)
pair <- map_chr(nominal_log, ~.x[[grep("Reading Expression files:",.x)]])
names(nominal_log) <- ss(ss(basename(pair),"_",2),"\\.")

time_lines <- map(nominal_log, ~.x[grep("time elapsed: ",.x)])
phenotype <- map_dbl(nominal_log, ~parse_number(.x[grep("\\^* \\d+ phenotypes",.x)][[1]]))
runtime <- map_dbl(nominal_log, ~max(parse_number(.x[grep("time elapsed: ",.x)])))

filter_log <- readLines(here("leafcutter","code", "logs", "08_filter_sQTL.txt"))

n_pairs <- parse_number(filter_log[grep("n pairs:", filter_log)])
n_pairs_FDR05 <- parse_number(gsub("n pairs FDR<0.05:","",filter_log[grep("n pairs FDR<0.05:", filter_log)]))

nominal_log_data <- tibble(feat ="Splice",
                           region = names(runtime),
                           n_feat = phenotype,
                           runtime,
                           n_pairs,
                           n_pairs_FDR05 )
# 
# # A tibble: 2 × 5
#   data     n_feat runtime   n_pairs n_pairs_FDR05
# <chr>     <dbl>   <dbl>     <dbl>         <dbl>
# 1 Amygdala 209476    19.2 451404371       7851177
# 2 sACC     209476    19.2 449487451       8324138

write.csv(nominal_log_data, file = here("leafcutter","data", "LC_sQTL_summary.csv"))
