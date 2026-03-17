## collect QDNAseq objects or data frames for single cell profiles into single object
# two examples usage cases
# Rscript merge.R target.rds input_folder file_in_folder
# Rscript merge.R target.rds input_1.rds input_2.rds input_3.rds

require(BiocGenerics, quietly=TRUE)
require(gtools, quietly=TRUE)
require(dplyr, quietly=TRUE)
#require(scAbsolute)
if (interactive()){
  BASEDIR="~/scAbsolute"
} else {
  BASEDIR="/opt/scAbsolute"
}
source(file.path(BASEDIR, "R/core.R"))
args = commandArgs(trailingOnly=TRUE)
indexSort = TRUE
options(future.globals.maxSize= 4096*1024^2)

rdsFiles = tail(args, n=-1)
if (length(args) >= 3){
  d = args[2]
  f = args[3]  
}else{
  d = args[2]
  f = NULL
}

filesToRemove <- character(0)

if (dir.exists(d) && (is.null(f) || !file.exists(f))){
  # merge based on filename ending (f) in folder (d) to file args[1]
  setwd(d)
  rdsFiles = list.files(pattern=paste0(f, ".rds$"), recursive=FALSE)
}else{
  # merge all files in args[2:end] to file args[1]
  stopifnot(all(endsWith(rdsFiles, ".rds")))
  stopifnot(all(file.exists(rdsFiles)))

  if(indexSort){
    rdsFiles = sort(rdsFiles, decreasing=FALSE)
  }

  # Remove truly empty files (process crashed before any output was written)
  for (file in rdsFiles) {
    size <- file.info(file)$size
    if (is.na(size) || size == 0) {
      file.remove(file)
      filesToRemove <- c(filesToRemove, file)
      cat("Empty/unreadable file removed (process crash):", file, "\n")
    }
  }
  rdsFiles <- setdiff(rdsFiles, filesToRemove)

}

if (is.data.frame(readRDS(rdsFiles[[1]]))) {
  rdsData = base::lapply(rdsFiles, readRDS)#, USE.NAMES = FALSE)
  mergedRDS = do.call(rbind, rdsData)
} else {
  rdsData = base::sapply(rdsFiles, readRDS, USE.NAMES = FALSE)  
  # alternative command, but issues with NaN values in pData: 
  #mergedRDS = do.call(BiocGenerics::combine, rdsData)
  mergedRDS = combineQDNASets(unlist(rdsData))
}

# Write failure report: crashed files + cells with failure_reason in pData
crashed_cells = data.frame(name=filesToRemove, failure_reason="process_crash", stringsAsFactors=FALSE)
if(!is.data.frame(mergedRDS) && "failure_reason" %in% colnames(Biobase::pData(mergedRDS))){
  pd = Biobase::pData(mergedRDS)
  pdata_failed = pd[!is.na(pd$failure_reason), c("name", "failure_reason")]
  failed_report = rbind(crashed_cells, pdata_failed)
}else{
  failed_report = crashed_cells
}
write.csv(failed_report, file=paste0(dirname(args[1]), "/failed_cells.csv"), row.names=FALSE)
saveRDS(mergedRDS, file = args[1])
