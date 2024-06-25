
library(here)
library(SingleCellExperiment)

# Dissection: Amygdaloid complex (AMY) - Central nuclear group - CEN
sce_amy <- Seurat::as.SingleCellExperiment(readRDS(here("deconvolution", "BICCN_data", "AMY_CEN", "c821776d-a9fd-4a9b-b0ef-65bee87c7dbc.rds")))
dim(sce_amy)
# 59236 42699

colData(sce_amy)
table(sce_amy$supercluster_term)
# Amygdala excitatory                           Astrocyte                     CGE interneuron 
#               3892                                2387                                4146 
# Committed oligodendrocyte precursor   Deep-layer corticothalamic and 6b       Deep-layer intratelencephalic 
#                                  42                                1144                                6981 
# Deep-layer near-projecting       Eccentric medium spiny neuron                           Ependymal 
#                         43                                3691                                  29 
# Fibroblast           LAMP5-LHX6 and Chandelier                     MGE interneuron 
#         51                                1386                                1760 
# Medium spiny neuron                           Microglia         Midbrain-derived inhibitory 
#                4708                                1951                                  20 
# Miscellaneous                     Oligodendrocyte           Oligodendrocyte precursor 
#           70                                5094                                3032 
# Splatter                 Thalamic excitatory                   Upper rhombic lip 
#     2092                                   1                                   9 
# Upper-layer intratelencephalic                            Vascular 
#                            109                                  61 

# Dissection: Cerebral cortex (Cx) - Subcallosal Gyrus (SCG) - Subgenual (subcallosal) division of MFC - A25
sce_sacc <- Seurat::as.SingleCellExperiment(readRDS(here("deconvolution", "BICCN_data", "SCG_MFC", "82faf671-658a-4c88-8a7d-618fc7f68fad.rds")))
table(sce_sacc$supercluster_term)
# Amygdala excitatory                           Astrocyte                     CGE interneuron 
#                 107                                 713                                4614 
# Committed oligodendrocyte precursor   Deep-layer corticothalamic and 6b       Deep-layer intratelencephalic 
#                                  23                                2005                                7204 
# Deep-layer near-projecting       Eccentric medium spiny neuron                          Fibroblast 
#                        909                                  12                                   9 
# Hippocampal CA1-3           LAMP5-LHX6 and Chandelier                     MGE interneuron 
#                 1                                1024                                6738 
# Medium spiny neuron                           Microglia                       Miscellaneous 
#                  32                                 361                                 385 
# Oligodendrocyte           Oligodendrocyte precursor                            Splatter 
#             804                                 732                                 105 
# Upper-layer intratelencephalic                            Vascular 
#                         11968                                  21

