library(QDNAseq)
args = commandArgs(trailingOnly=TRUE)

BASEDIR= args[1]

source(file.path(BASEDIR, "R/core.R"))
source(file.path(BASEDIR, "R/visualization.R"))

CN = readRDS(args[2])

# plot single cell copy number profile
plotCopynumber(CN[, 1])

# plot chromosome
chr1 = startsWith(x = rownames(CN), prefix="1:")
plotCopynumber(CN[chr1, 1])

# plot only CN profile (e.g. for publication)
chr1 = startsWith(x = rownames(CN), prefix="1:")
plotCopynumber(CN[, 1],showUnique = FALSE, main = "", readinfo = FALSE, showMarker = FALSE)
