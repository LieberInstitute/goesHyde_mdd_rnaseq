library("tidyverse")
library("jaffelab")
library("ggrepel")
library("here")

nominal_log_fn <- list.files(here("eqtl","code", "logs"), pattern = "tensorqtl_genomewide_nominal_gpu*", full.names = TRUE)

nominal_log <- map(nominal_log_fn, readLines)
pair <- map_chr(nominal_log, ~.x[[grep("Reading Expression files:",.x)]])
names(nominal_log) <- ss(ss(pair,"/", 5),"\\.")

time_lines <- map(nominal_log, ~.x[grep("time elapsed: ",.x)])
phenotype <- map_dbl(nominal_log, ~parse_number(.x[grep("\\^* \\d+ phenotypes",.x)][[1]]))
runtime <- map_dbl(nominal_log, ~max(parse_number(.x[grep("time elapsed: ",.x)])))

filter_log <- readLines(here("eqtl","code", "logs", "filter_genomewide_nominal.txt"))
n_pairs <- parse_number(filter_log[grep("n pairs:", filter_log)])
n_pairs_FDR05 <- parse_number(gsub("n pairs FDR<0.05:","",filter_log[grep("n pairs FDR<0.05:", filter_log)]))

pairs_v_time <- tibble(data = names(runtime),
                       n_feat = phenotype,
                       runtime,
                       n_pairs,
                       n_pairs_FDR05 )

# data          n_feat runtime   n_pairs n_pairs_FDR05
# <chr>          <dbl>   <dbl>     <dbl>         <dbl>
# 1 gene_Amygdala  25132    2.52  53637238        771221 * 25212
# 2 gene_sACC      25132    2.55  53403152        853567
# 3 exon_Amygdala 399023   45.6  840508027       5960604 * 399529
# 4 exon_sACC     399023   45.6  836725998       7192434
# 5 jxn_Amygdala  304125   37.0  633281109       9381590 * 305558
# 6 jxn_sACC      304125   36.9  630534075      10434226
# 7 tx_Amygdala    79565    8.21 167717620       1684554 * 79730
# 8 tx_sACC        79565    8.22 166980983       1906150

pairs_v_time_scatter <- ggplot(pairs_v_time, aes(x = n_pairs, y = runtime))+
  geom_point() +
  geom_text_repel(aes(label = data)) +
  labs(y = "Runtime (min)")

ggsave(pairs_v_time_scatter, filename = here("eqtl", "plots", "test", "pairs_v_time_scatter.png"))

## checks 
pairs_v_time %>%
  group_by(data) %>%
  transmute(prop = n_pairs_FDR05/n_pairs)

pairs_v_time %>%
  group_by(data) %>%
  transmute(mean_snp = n_pairs/n_feat)

#### Significant barplot ####
n_signif <- pairs_v_time %>%
  mutate(n_pairs_ns = n_pairs - n_pairs_FDR05) %>%
  select(data, n_pairs_ns, n_pairs_FDR05) %>%
  pivot_longer(!data, names_to = "Signif", values_to = "n", names_prefix ="n_pairs_")

signif_barplot <- n_signif %>%
  ggplot(aes(x = data, y = n, fill = Signif)) +
  geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(signif_barplot, filename = here("eqtl", "plots", "test", "signif_barplot.png"))


