## run script for scUnique
# scUnique - pipeline with joint copy number inference, FDR correction and rCNA calling for scDNAseq data

## Variables ####
args = commandArgs(trailingOnly=TRUE)
set.seed(2020)

if (interactive()){
  rm(list=ls())
  
  # it is possible to combine multiple files like this
  sampleFile = c(
    "~/Data/download/scDNAseq-workflow/results/500/PEO1_500.rds",
    "~/Data/download/scDNAseq-workflow/results/500/PEO4_500.rds"
  )
  
  sampleName = "PEO1-PEO4"
  RESULTPATH = paste0("~/Data/download/scDNAseq-workflow/results/", "")
  
  suppressPackageStartupMessages({
    require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)
    require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
  })
  
  reticulate::use_condaenv(condaenv = "conda_runtime", conda = "/opt/conda/bin/conda")
  .libPaths("/opt/conda/envs/conda_runtime/lib/R/library")
  BASEDIR="~/"
  breakEarly=FALSE
  change_prob = NULL
  
  # sample subset of data for debugging purposes
  sampleFile = do.call("c", (base::strsplit(sampleFile, split=",")))
  
  if(length(sampleFile) == 1){
    CN = readRDS(sampleFile) 
  }else{
    source(file.path(BASEDIR, "scAbsolute/R/core.R"))
    CN = combineQDNASets(future.apply::future_lapply(sampleFile, readRDS))
  }
  
}else{
  
  sampleFile = args[1]
  RESULTPATH = args[2]
  change_prob = as.numeric(args[3])
  
  suppressPackageStartupMessages({
    require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)
    require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
  })
  reticulate::use_condaenv(condaenv = "conda_runtime", conda = "/opt/conda/bin/conda")
  BASEDIR="/opt"
  
  # sample subset of data for debugging purposes
  sampleFile = do.call("c", (base::strsplit(sampleFile, split=",")))
  
  if(length(sampleFile) == 1){
    CN = readRDS(sampleFile) 
  }else{
    source(file.path(BASEDIR, "scAbsolute/R/core.R"))
    CN = combineQDNASets(future.apply::future_lapply(sampleFile, readRDS))
  }
}

BASEDIR=normalizePath(BASEDIR)
require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)
Sys.setenv(RETICULATE_PYTHON="/opt/conda/envs/conda_runtime/bin/python")
Sys.setenv(TENSORFLOW_PYTHON="/opt/conda/envs/conda_runtime/bin/python")
reticulate::use_condaenv("conda_runtime")
# segmentation.py is in scAbsolute
#reticulate::source_python(file.path("/opt/scAbsolute/R/segmentation.py"), convert=TRUE)
require(tensorflow)

# PACKAGE DEPENDENCIES ====
suppressPackageStartupMessages({
  
  require(devtools, quietly=TRUE, warn.conflicts = FALSE)
  require(future.apply, quietly=TRUE, warn.conflicts = FALSE)
  plan(multicore)
  print(paste0("available cores: ", availableCores()))
  require(progressr, quietly=TRUE, warn.conflicts = FALSE)
  handlers("progress")
  require(fastcluster, quietly=TRUE, warn.conflicts = FALSE)
  require(phylogram, quietly=TRUE, warn.conflicts = FALSE)
  require(dendextend, quietly=TRUE, warn.conflicts = FALSE)
  require(tiff, quietly=TRUE, warn.conflicts = FALSE)
  require(stats, quietly=TRUE, warn.conflicts = FALSE)
  require(dynamicTreeCut, quietly=TRUE, warn.conflicts = FALSE)
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
  
  require(Biobase, quietly = TRUE, warn.conflicts = FALSE)
  require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
  require(cowplot, quietly = TRUE, warn.conflicts = FALSE)
  require(matrixStats, quietly = TRUE, warn.conflicts = FALSE)
  
  require(parallel, quietly = TRUE, warn.conflicts = FALSE)
  require(foreach, quietly = TRUE, warn.conflicts = FALSE)
  require(doParallel, quietly = TRUE, warn.conflicts = FALSE)
  options(future.globals.onReference = "error")
  
  # scAbsolute dependencies
  source(file.path(BASEDIR, "scAbsolute/R/core.R"))
  # import combineQDNASets, binsToUseInternal
  source(file.path(BASEDIR, "scAbsolute/R/visualization.R"))
  # visualize individual cells
  
  source(file.path(BASEDIR, "scUnique/R/core.R"))
  source(file.path(BASEDIR, "scUnique/R/segmentation.R"))
  source(file.path(BASEDIR, "scUnique/R/scUnique.R"))
  source(file.path(BASEDIR, "scUnique/R/visualize.R"))
  source(file.path(BASEDIR, "scUnique/R/postprocessing.R"))
  
  use_condaenv("conda_runtime")
})
## Detect non-exportable objects and give an error asap
# PACKAGE DEPENDENCIES END
ncores <- Sys.getenv('NCORES')
if(ncores == ""){ncores = 8}
ncores = max(as.integer(ncores)-2, 4)
print(paste0("Using ", ncores, " cores."))
stopifnot(ncores >= 4)

## PARAMETERS ====
randomSeed = 2020
max_state = 8
pvalue_cutoff = 0.01
n_samples_lrtest=10000
start_time_all <- Sys.time()

# size of events to consider (inclusive)
minimum_eventsize=3
# size of absolute difference in copynumber to override minimum eventsize (inclusive)
focalDifference=2
if(is.null(change_prob)){
  change_prob = 0.05
}
# if both set to 0, do not filter medicc segments
mediccSizeCutoff=0#5
mediccDiffAmplitude=0#2

print(paste0("Change prob: ", change_prob))
windowSize=0
## END PARAMETERS

df = Biobase::pData(CN) %>%  tidyr::separate(name, sep="_", into=c("UID", "SLX", "cellid", "celltag"), remove=FALSE)
print(paste0("SAMPLE FILES: ", sampleFile))
stopifnot(all(df$splitPerChromosome) || all(!df$splitPerChromosome))
splitPerChromosome=df$splitPerChromosome[1]


# QC filtering is done prior to this analysis and cells passing quality inspection are described in 
pass_qc = do.call("c", lapply(list.files(path=file.path(RESULTPATH, "pass_qc/"), pattern="\\.tsv", full.names = TRUE),
                  function(x){readr::read_tsv(x, col_names=c("UID", "SLX", "name"), col_types="ccc")[["name"]]}))
include_cells = colnames(CN)[colnames(CN) %in% pass_qc]
if(length(include_cells)==0){
  warning("It seems like no cell passed qc.\nDid you create a pass_qc file for the sample?")
  stop()
}
# we generally don't trust Y and X chromosome calls for rCNA analysis
include_chr = !(startsWith(rownames(CN), "Y:"))
object = CN[include_chr, include_cells]
print(paste0("Total number of cells: ", dim(object)[2]))
print(paste0("Total number of included cells: ", length(include_cells), " (", format(round(length(include_cells)/(dim(CN)[2]), 3), nsmall = 2), ")"))
stopifnot(all((Biobase::pData(object)[["rpc"]] > 0 & Biobase::pData(object)[["hmm.alpha"]] > 0)))


## scUnique calls ====
df = Biobase::pData(object)
stopifnot(all(df$rpc > 0))
valid=binsToUseInternal(object)
stopifnot(all(Biobase::assayDataElement(object, "copynumber")[valid,] <= max_state))

## joint breakpoint inference - using MEDICC2 ----
n_cells = dim(object)[[2]]
# NOTE this should be excluding
HMMPATH = dirname(dirname(Biobase::lcPrefix(sampleFile)))
df = getMediccTable(object, mediccSizeCutoff=mediccSizeCutoff, mediccDiffAmplitude=mediccDiffAmplitude)
df_alt = getMediccTable(object, mediccSizeCutoff=0, mediccDiffAmplitude=0)

## Refine segmentation ====
start_time_jointInference <- Sys.time()
jointSegmentationObject = jointSegmentationUpdateMEDICC(object,
                                                        change_prob = change_prob, windowSize=windowSize,
                                                        mediccSizeCutoff=mediccSizeCutoff, mediccDiffAmplitude=mediccDiffAmplitude,
                                                        splitPerChromosome=splitPerChromosome, ncores=ncores, HMMPATH=HMMPATH, BASEDIR)
end_time_jointInference <- Sys.time()
print(paste0("Run joint segmentation algorithm ",  difftime(end_time_jointInference,start_time_jointInference,units="mins")))

newObject = jointSegmentationObject$newObject
pdms_initial = jointSegmentationObject$pdms
tree_initial = jointSegmentationObject$medicc_tree


## MEDICC analysis - on fully segmented data ====
df = getMediccTable(newObject, mediccSizeCutoff=mediccSizeCutoff, mediccDiffAmplitude=mediccDiffAmplitude) # NOTE, by default accessing copynumber slot!
print("Segments for MEDICC run")
print(nrow(df %>% dplyr::filter(sample_id == unique(df$sample_id)[[1]])))
print("Samples for MEDICC run")
print(dim(newObject)[[2]])

require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
reticulate::py_discover_config()
e = new.env()
reticulate::source_python(file.path(BASEDIR, "scUnique/interface_medicc2.py"), convert=TRUE, envir=e)

start_time_medicc <- Sys.time()
out = reticulate::py_capture_output({
results = e$interface_medicc(df, as.integer(ceiling(ncores * (2/3))))
})
cat(out)
end_time_medicc <- Sys.time()
print(paste0("Run MEDICC algorithm 1 ",  difftime(end_time_medicc,start_time_medicc,units="mins")))

tree = results[[1]]; input_df = results[[2]]; output_df = results[[3]];
summary = results[[4]]; pdms = results[[5]];
# unique events is difference between leaf node and MRA (first node up in the tree)
uniqueEvents = results[[6]]; newick_tree = results[[7]];

medicc_tree <- read.dendrogram(text = newick_tree)
# NOTE remove diploid
tree_pre <- phylogram::prune(medicc_tree, pattern = "diploid")

## Create background profile
valid = binsToUseInternal(newObject)
dims = dim(newObject@assayData$copynumber[valid,])

unique.mat = matrix(0, nrow=dims[[1]], ncol=dims[[2]])
profiles_background = newObject@assayData$copynumber
profiles_background_initial = newObject@assayData$copynumber
colnames(unique.mat) = colnames(newObject)
colnames(profiles_background) = colnames(newObject)
cellnames = colnames(newObject)

counter = 0
for(cellname in cellnames){
  # print(cellname)
  ue = uniqueEvents[[cellname]][["events"]]
  if(length(ue) == 0){
    next
  }
  for(uindex in 1:length(ue)){
    # print(uindex)
    event = ue[[uindex]]
    if(length(event) == 0){
      next
    }
    
    counter = counter + event$end - event$start + 1
    
    unique.mat[event$start:event$end, cellname] = 1
    # stopifnot(all(min(profiles_background[valid,cellname][event$start:event$end], max_state) == as.numeric(event$h1)))
    if(event$chr == "23"){
      event$chr = "X"
    }
    if(event$chr == "24"){
      event$chr = "Y"
    }
    # stopifnot(all(base::startsWith(names(profiles_background[valid,cellname][event$start:event$end]), event$chr)))
    
    profiles_background[valid,cellname][event$start:event$end] = min(as.numeric(event$h0), max_state)
    counter = counter + (event$end - event$start)
  }
  
}

start_time_probeEvidence = Sys.time()
print(paste0("Run probe evidence (", as.POSIXlt(start_time_probeEvidence), ")"))
evidence = probe_evidence(newObject, profiles_background,
                          pvalue_cutoff=pvalue_cutoff, minimum_eventsize=minimum_eventsize, focalDifference = focalDifference,
                          n_samples=n_samples_lrtest, splitPerChromosome=splitPerChromosome, ncores=ncores, HMMPATH=HMMPATH, BASEDIR=BASEDIR)

end_time_probeEvidence = Sys.time()
print(paste0("Completed probe evidence (", as.POSIXlt(end_time_probeEvidence), ")"))


newObjectControl = evidence$object
df_events = evidence$df_events
df_pass = evidence$df_pass
df_rejected = evidence$df_rejected
# objects using only low frequency rCNAs
newObjectControlFreq = evidence$objectFreq
profiles_background_freq = evidence$backgroundFreq


df = getMediccTable(newObjectControl, mediccSizeCutoff=0, mediccDiffAmplitude=0) # NOTE, by default accessing copynumber slot!
print("Segments for MEDICC run")
print(nrow(df %>% dplyr::filter(sample_id == unique(df$sample_id)[[1]])))
print("Samples for MEDICC run")
print(dim(newObject)[[2]])

require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
reticulate::py_discover_config()
reticulate::source_python(file.path(BASEDIR, "scUnique/interface_medicc2.py"), convert=TRUE)

start_time_medicc <- Sys.time()
out = reticulate::py_capture_output({
results = interface_medicc(df, as.integer(ceiling(ncores * (2/3))))
})
cat(out)
end_time_medicc <- Sys.time()
print(paste0("Run MEDICC algorithm 2 ",  difftime(end_time_medicc,start_time_medicc,units="mins")))

tree = results[[1]]; input_df = results[[2]]; output_df = results[[3]];
summary = results[[4]]; pdms = results[[5]];
# unique events is difference between leaf node and MRA (first node up in the tree)
uniqueEvents = results[[6]]; newick_tree = results[[7]];

# format tree_profiles: list(name, ancestor.name, seq.short, seq.long, isLeafNode)
tree_profiles = results[[8]];

medicc_tree <- read.dendrogram(text = newick_tree)
# NOTE remove diploid
tree_post <- phylogram::prune(medicc_tree, pattern = "diploid")

# compute df_rcna
df_rcna = simple_rcnas(pdms, newObjectControl, ncores=ncores)

# postprocessing
df_rcna_post = metainformation_rCNA(newObjectControl, df_rcna, BASEDIR=BASEDIR) 
df_pass_post = metainformation_rCNA(newObjectControl, df_pass, BASEDIR=BASEDIR) 

## number of events
n_events = dplyr::right_join(df_pass %>% dplyr::count(cellname, .drop=FALSE) %>% dplyr::ungroup(), dplyr::tibble(cellname=colnames(object)), by="cellname")
n_events$n = ifelse(is.na(n_events$n), 0, n_events$n)

## List of things to export
# medicc tree, post tree, different objects, pre tree
# df_pass, n_events, df_events
# profile_background

# newObjectControl == finalCN
# NOTE, final profiles only keeping TRUE events -> profiles_unique == copynumber slot
# segmented slot: background profile (H0)
# copynumber slot: H0 + validated unique events
print("FINAL CN COMPLETED")
end_time_all <- Sys.time()


##########################################################
## Save data for later processing
FILENAME = paste0(lapply(sampleFile, function(x) sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(x))), collapse = "+")
if(!dir.exists(dirname(file.path(RESULTPATH, FILENAME, "a")))) dir.create(dirname(file.path(RESULTPATH, FILENAME, "a")), recursive = TRUE)
## Format RESULTPATH with sample ID
saveRDS(CN, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".rawCN.RDS")))
saveRDS(object, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".inputCN.RDS")))
saveRDS(newObject, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".jointSegmentedCN.RDS")))
saveRDS(newObjectControl, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".finalCN.RDS")))
saveRDS(newObjectControlFreq, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".finalCNfreq.RDS")))

saveRDS(tree_pre, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".tree_pre.RDS")))
saveRDS(medicc_tree, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".medicc_tree.RDS")))
saveRDS(tree_post, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".tree.RDS")))
saveRDS(uniqueEvents, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".uniqueEvents.RDS")))
saveRDS(pdms, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".pdms")))
saveRDS(phylogram::as.dendrogram(tree_post), file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".dendrogram.RDS")))
saveRDS(tree_profiles, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".tree_profiles.RDS")))

saveRDS(profiles_background, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".profiles_background.RDS")))
saveRDS(profiles_background_initial, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".profiles_background_initial.RDS")))
saveRDS(profiles_background_freq, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".profiles_background_freq.RDS")))

saveRDS(df_events, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".df_events.RDS")))
saveRDS(df_pass, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".df_pass.RDS")))
saveRDS(df_rcna, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".df_rcna.RDS")))
saveRDS(df_pass_post, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".df_pass_post.RDS")))
saveRDS(df_rcna_post, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".df_rcna_post.RDS")))
saveRDS(df_rejected, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".df_rejected.RDS")))
saveRDS(n_events, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".n_events.RDS")))
saveRDS(results, file=file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".medicc_results.RDS")))

save.image(file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".RData")))
##########################################################

print("DEBUGGING INFO:")
## Reproducibility info
devtools::session_info()
print("------------------------------")

## WRITE OUTPUT TO DEBUG FILE
capture.output({
  print("------------------------------")
  print("PROCESSING INFORMATION:")
  print("------------------------------")
  ## Some run information
  print(paste0("Runtime MEDICC algorithm ",  difftime(end_time_medicc,start_time_medicc,units="mins")))
  print(paste0("Runtime MEDICC joint algorithm ",  difftime(end_time_jointInference,start_time_jointInference,units="mins")))
  print(paste0("Runtime probe evidence ",  difftime(end_time_probeEvidence,start_time_probeEvidence,units="mins")))
  print(paste0("Runtime program ",  difftime(end_time_all,start_time_all,units="mins")))
  print("-------")
  # Copy number state distribution
  tbl = table(newObjectControl@assayData$copynumber)
  print("Distribution of copy number states across sample")
  print(prop.table(tbl))
  print("-------")
  # Issues with sample:
  print("Total number of cells kept")
  print(paste0(dim(newObjectControl)[[2]], " cells (out of ", dim(CN)[[2]], " cells)"))
  print(paste0("Percentage of cells discarded: ", (dim(CN)[[2]]-dim(newObjectControl)[[2]]) / dim(CN)[[2]]))
  print("-------")
  print("-------")
  print(sampleFile)
  print("END")},
  file = file.path(RESULTPATH, FILENAME, paste0(FILENAME, ".txt")))

print("COMPLETED")
