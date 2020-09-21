library(sessioninfo)

### lets start w/ hg19 for the remaining ~3M
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
message("\n ***** hg19 coordinates ***** ", Sys.time())
dbSnp142_list = mclapply(paste0("ch",c(1:22,"X")), function(x) {
  cat(x,", ")
  y = SNPlocs.Hsapiens.dbSNP142.GRCh37::getSNPlocs(x, as.GRanges=TRUE)
  seqlevels(y) = gsub("ch", "chr", seqlevels(y))
  y$RefSNP_id = paste0("rs", y$RefSNP_id)
  return(y)
},mc.cores=6)
dbSnp142 = unlist(GRangesList(dbSnp142_list))

message("dbSnp142 lenght: ",length(dbSnp142))

############# hg38 coordinates ######################
message("\n ***** hg38 coordinates ***** ", Sys.time())
library(SNPlocs.Hsapiens.dbSNP149.GRCh38)
dbSnp149_list = mclapply(c(1:22,"X"), function(x) {
  cat(x,", ")
  y = snpsBySeqname(SNPlocs.Hsapiens.dbSNP149.GRCh38, x)
  seqlevels(y) = paste0("chr", seqlevels(y))
  return(y)
},mc.cores=6)

dbSnp149 = unlist(GRangesList(lapply(dbSnp149_list, as, "GRanges")))
message("dbSnp149 lenght: ",length(dbSnp149))

save(dbSnp142, dbSnp149, file = "dbSNP.Rdata")

# sgejobs::job_single('get_dbSNP', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript get_dbSNP.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
