library(QDNAseq)
args = commandArgs(trailingOnly=TRUE)

BASEDIR= args[1] #/opt/scAbsolute #or #/rds/user/adr44/hpc-work/scAbsolute/github/scAbsolute/

source(file.path(BASEDIR, "R/core.R"))
source(file.path(BASEDIR, "R/visualization.R"))

CN = readRDS(args[2]) #example: #/rds/user/adr44/hpc-work/scAbsolute/github/scDNAseq-workflow/mouse_compatibility/results/scale/500/individual/SH171012_I_028.rds

# plot single cell copy number profile
plotCopynumber(CN[, 1])

# plot chromosome
chr1 = startsWith(x = rownames(CN), prefix="1:")
plotCopynumber(CN[chr1, 1])

# plot only CN profile (e.g. for publication)
chr1 = startsWith(x = rownames(CN), prefix="1:")
plotCopynumber(CN[, 1],showUnique = FALSE, main = "", readinfo = FALSE, showMarker = FALSE)
dev.off()