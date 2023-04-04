library(QDNAseq)
BASEDIR="~/Projects/scAbsoluteMatters/scAbsolute"
source(file.path(BASEDIR, "R/core.R"))
source(file.path(BASEDIR, "R/visualization.R"))
CN = readRDS("~/Projects/scAbsoluteMatters/results/500/individual/UID-10X-Fibroblast-cell_SLX-00000_000001_AAACCTGAGCAGCCTC-1.rds")

# plot single cell copy number profile
plotCopynumber(CN[, 1])

# plot chromosome
chr1 = startsWith(x = rownames(CN), prefix="1:")
plotCopynumber(CN[chr1, 1])

# plot only CN profile (e.g. for publication)
chr1 = startsWith(x = rownames(CN), prefix="1:")
plotCopynumber(CN[, 1],showUnique = FALSE, main = "", readinfo = FALSE, showMarker = FALSE)