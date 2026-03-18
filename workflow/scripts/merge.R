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

  # Repair any per-cell objects where assayData colnames diverged from pData
  # rownames (can happen silently via pData<- assignments in scAbsolute).
  # Use pData$name as the canonical cell identifier set by QDNAseq.
  rdsData = lapply(rdsData, function(obj) {
    if (!is(obj, "QDNAseqCopyNumbers")) return(obj)
    canonical <- Biobase::pData(obj)[["name"]]
    adata_names <- colnames(Biobase::assayDataElement(obj, "copynumber"))
    pdata_rownames <- rownames(Biobase::pData(obj))
    if (!identical(adata_names, pdata_rownames)) {
      cat("Repairing sampleName mismatch for cell:", canonical, "\n")
      cat("  assayData colname :", adata_names, "\n")
      cat("  pData rowname     :", pdata_rownames, "\n")
      # Fix pData rownames to match assayData colnames (QDNAseq sets both
      # to basename(bam) without .bam; assayData colnames are more reliable).
      pd <- Biobase::pData(obj)
      rownames(pd) <- adata_names
      Biobase::pData(obj) <- pd
    }
    return(obj)
  })

  # alternative command, but issues with NaN values in pData:
  #mergedRDS = do.call(BiocGenerics::combine, rdsData)
  mergedRDS = tryCatch(
    combineQDNASets(rdsData),
    error = function(e) {
      if (grepl("sampleNames differ", conditionMessage(e))) {
        # Replicate what combineQDNASets does to find the mismatch
        qdna_list <- Filter(function(o) is(o, "QDNAseqCopyNumbers"), rdsData)
        combined_adata_names <- unlist(lapply(qdna_list, function(o)
          colnames(Biobase::assayDataElement(o, "copynumber"))))
        combined_pdata_names <- unlist(lapply(qdna_list, function(o)
          rownames(Biobase::pData(o))))
        cat("combineQDNASets failed: sampleNames differ\n")
        cat("  combined assayData colnames (first 5):",
            paste(head(combined_adata_names, 5), collapse=", "), "\n")
        cat("  combined pData rownames   (first 5):",
            paste(head(combined_pdata_names, 5), collapse=", "), "\n")
        mismatch_idx <- which(!mapply(identical, combined_adata_names, combined_pdata_names))
        if (length(mismatch_idx) > 0) {
          cat("  mismatched positions:", paste(head(mismatch_idx, 10), collapse=", "), "\n")
          for (i in head(mismatch_idx, 5)) {
            cat(sprintf("    [%d] assayData='%s'  pData rowname='%s'\n",
                        i, combined_adata_names[i], combined_pdata_names[i]))
          }
        } else {
          # names match per-element but identical() on full vectors fails — check lengths
          cat("  lengths: assayData", length(combined_adata_names),
              " pData", length(combined_pdata_names), "\n")
          cat("  pData ncol per cell (first 5):",
              paste(head(sapply(qdna_list, function(o) ncol(Biobase::pData(o))), 5), collapse=", "), "\n")
        }
      }
      stop(e)
    }
  )
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
