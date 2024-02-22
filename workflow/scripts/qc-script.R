## script for interactive QC analysis of scAbsolute output
# We recommend manually running and inspecting the results,
# while it's theoretically possible to run this script automated as well

# NOTE: by default, datasets are split by UID and SLX variable,
# so it is important to set these two variables in the df data frame
# we throw an error below if this hasn't happened

## Variables ####
args = commandArgs(trailingOnly=TRUE)
set.seed(2020)

if (interactive()){
  rm(list=ls())
  
  sampleFile = c(
    "~/scDNAseq-workflow/results/500/PEO1_500.rds",
    "~/scDNAseq-workflow/results/500/PEO4_500.rds"
    
  )
  sampleName = "PEO1-PEO4"
  
}else{
  sampleName = args[1]
  sampleFile = args[2]
}
sampleFile = do.call("c", (base::strsplit(sampleFile, split=",")))

## Setup ====
BASEDIR="~/"
BASEDIR=normalizePath(BASEDIR)
WORKFLOW_PATH="~/scDNAseq-workflow/"
require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)
require(ggbeeswarm, quietly = TRUE, warn.conflicts = FALSE)
require(ggpubr, quietly = TRUE, warn.conflicts = FALSE)
require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
require(ComplexHeatmap, quietly = TRUE, warn.conflicts = FALSE)
source(file.path(BASEDIR, "scUnique/R/core.R"))
source(file.path(BASEDIR, "scUnique/R/visualize.R"))
source(file.path(BASEDIR, "scAbsolute/R/core.R"))
source(file.path(BASEDIR, "scAbsolute/R/visualization.R"))

object = combineQDNASets(future.apply::future_lapply(sampleFile, readRDS))

# Check for basic properties - here we use the common name format to directly create UID/SLX from filenames
df = Biobase::pData(object) %>%
 tidyr::separate(name, sep="_", into=c("UID", "SLX", "cellid", "celltag"), remove=FALSE)
stopifnot("UID" %in% colnames(df))
stopifnot("SLX" %in% colnames(df))
 
# 1. number of reads extreme outliers
reads_cutoff = quantile(df$used.reads, c(0.025, 0.975)) # by default we remove the top and bottom 2.5% of cells with very high or low read number 
p_coverage = ggplot(data=df %>% dplyr::select(rpc, used.reads, UID, SLX)) +
  geom_histogram(aes(x=used.reads, y=..density.., color=interaction(UID, SLX)), bins = 100, position="identity", fill="white") +
  geom_density(aes(x=used.reads, color=interaction(UID, SLX))) +
  geom_vline(xintercept = reads_cutoff, color="red") +
  theme_pubclean()
issue_coverage = colnames(object)[(Biobase::pData(object)[["used.reads"]] < reads_cutoff[[1]]) | 
                                  (Biobase::pData(object)[["used.reads"]] > reads_cutoff[[2]])]
p_coverage
if(any(Biobase::pData(object)[["rpc"]] < 25)){
  warning("Some cells have an rpc values of < 25!\nThis might indicate that you are under powered to detect rCNAs and copy number alterations at this bin size.")
  print(Biobase::pData(object)[["name"]][Biobase::pData(object)[["rpc"]] < 25]) 
}

# manually add cells that should be excluded
issue_manual = c()
exclude_cells = base::union(issue_coverage, issue_manual)
stopifnot(all(is.character(exclude_cells)))
include = setdiff(colnames(object), exclude_cells)

object = object[, include]
df = df %>% dplyr::filter(!(name %in% base::union(issue_coverage, issue_manual)))
stopifnot(dim(df)[[1]] == dim(object)[[2]])
print(paste0("ISSUES: ", base::union(issue_coverage, issue_manual)))

# PARAMETERS ====
# see source for description
# replicating
cutoff_replicating = 1.50
cutoff_replicating_iqr = 1.5
cutoff_replicating_hard = NULL
# MAPD
mapd_cutoff = 2.0
mapd_density_control = FALSE
# Gini
gini_norm_cutoff = 2.00
gini_density_control = FALSE
# alpha
alpha_cutoff = 1.5
alpha_hard_cutoff = 0.05 # 0.05: about 99 percentile of data for 500 kb
print(paste0("Number of cells: ", dim(object)[[2]]))

## Identify replicating and outlier cells ====


# 2. S phase cells
df = predict_replicating(df %>%
 dplyr::mutate(UID2 = UID, SLX2=SLX),
 cutoff_value=cutoff_replicating, iqr_value = cutoff_replicating_iqr, hard_threshold = cutoff_replicating_hard)
fig_cellcycle = ggplot(data = df) +
  geom_quasirandom(aes(x=SLX, y=cycling_activity, color=replicating, name=name), size=0.8, alpha=1.0) +
  geom_hline(aes(yintercept=cycling_cutoff_sd)) +
  geom_hline(aes(yintercept=cycling_median), linetype="dashed") +
  facet_wrap(~UID, scales = "free_x") + theme_pubclean()
fig_cellcycle

issue_sphase = df$name[df$replicating]

df$cycling = df$name %in% issue_sphase
table(df$cycling)
prop.table(table(df$cycling, df$UID), margin = 2)


# 3. QC control (after removal of cycling cells)
iq = qc_alpha(
              qc_gini_norm(
                            qc_mapd(df %>% dplyr::filter(!replicating), 
                            cutoff_value = mapd_cutoff, densityControl=mapd_density_control),
              cutoff_value = gini_norm_cutoff, densityControl = gini_density_control),
    cutoff_value=alpha_cutoff, hard_cutoff = alpha_hard_cutoff)
iq$outlier = iq$dmapd.outlier | iq$dgini.outlier | iq$alpha.outlier
aq = dplyr::left_join(df, iq, by="name", suffix=c("", ".y")) %>% dplyr::mutate(outlier = case_when(
                                              replicating ~ "replicating",
                                              name %in% iq$name[iq$outlier] ~ "outlier",
                                              TRUE ~ "pass"))
  
alpha_value = 0.7
scale_values = c(alpha_value, 0.1, alpha_value, alpha_value)
q1 = ggplot(data = aq) + geom_point(aes(x=rpc, y=dmapd, color=outlier, name=name, alpha=outlier)) + facet_wrap(~UID+SLX, scales = "free_y", nrow=1) + theme_pubclean() + scale_alpha_manual(values = scale_values)
q2 = ggplot(data = aq) + geom_point(aes(x=rpc, y=dgini, color=outlier, name=name, alpha=outlier)) + facet_wrap(~UID+SLX, scales = "free_y", nrow=1) + theme_pubclean() + scale_alpha_manual(values = scale_values)
q4 = ggplot(data = aq) + geom_point(aes(x=rpc, y=n_segs, color=outlier, name=name, alpha=outlier)) + facet_wrap(~UID+SLX, scales = "free_y", nrow=1) + theme_pubclean() + scale_alpha_manual(values = scale_values) 
q5 = ggplot(data = aq) + geom_quasirandom(aes(x=NA, y=error_seg_l2, color=outlier, name=name, alpha=outlier)) + facet_wrap(~UID+SLX, scales = "free_y", nrow=1) + theme_pubclean() + scale_alpha_manual(values = scale_values)
p3 = ggplot(data = aq %>% dplyr::filter(!replicating, rpc > 15)) + geom_quasirandom(aes(x=alpha.outlier, y=hmm.alpha, color=dmapd.outlier | dgini.outlier, name=name, alpha=outlier)) + geom_hline(yintercept=alpha_hard_cutoff, color="red", linetype="dashed") + facet_wrap(~UID+SLX, nrow=1) + theme_pubclean() + scale_alpha_manual(values = c(alpha_value, 0.1, alpha_value, 0.9, 1.0, 1.0)) + scale_y_continuous(trans="log10")

fig_qc = ggpubr::ggarrange(
                  ggpubr::ggarrange(q1  + theme(legend.position = "none") + rremove("xlab"), #, strip.background = element_blank(), strip.text = element_blank()),
                                    q2  + theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) + rremove("xlab"),
                                    q4  + theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) + rremove("xlab"),
                                    q5 + theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank()) + rremove("xlab") + 
                                      guides(colour = guide_legend(nrow=1, override.aes = list(alpha = 1, size=2))) + 
                                      theme(legend.key = element_rect(fill = "white"),
                                            legend.title = element_text(size = 8, face = 2))
                                      , ncol=2, nrow=2),
                  p3,
                  ncol=1, heights = c(2, 1))
fig_qc 
issue_qc = iq$name[iq$outlier]
df$remove = df$name %in% c(issue_qc)
table(df$remove)
prop.table(table(df$remove, df$UID), margin = 2)


# 4. Ploidy (cells with known wrong ploidy can be removed)
# consider re-running scAbsolute with better or more specific
# ploidy constraints if at all possible
p_ploidy = ggplot(data = aq) +
  geom_density2d(aes(x=rpc, y=ploidy), color="black") + 
  geom_point(aes(x=rpc, y=ploidy, color=outlier, name=name)) + theme_pubclean() + scale_alpha_manual(values = scale_values)
p_ploidy2 = ggplot(data = aq %>% dplyr::filter(rpc > 25)) + geom_density(aes(x=ploidy, color=outlier)) + theme_pubclean() # + facet_wrap(~UID)
p_ploidy2_low = ggplot(data = aq %>% dplyr::filter(rpc < 25)) + geom_density(aes(x=ploidy, color=outlier)) + theme_pubclean() 
p_reads = ggplot(data = aq) + geom_quasirandom(aes(x=interaction(UID, SLX), y=rpc, color=outlier, name=name, alpha=outlier)) + theme_pubclean() + scale_alpha_manual(values = scale_values) + geom_hline(yintercept=25, color="red", linetype="dashed") + theme(axis.text.x = element_text(angle=270))
fig_ploidy = ggpubr::ggarrange(p_ploidy, ggpubr::ggarrange(p_ploidy2, p_ploidy2_low, nrow=2), p_reads + theme(legend.position = "none"), nrow=1)
fig_ploidy

issue_ploidy = base::union(iq$name[iq$ploidy < 1.1], iq$name[iq$ploidy > 8.0])


## Optionally visualize copy number profiles
names = colnames(object)
pass_names = names[!(names %in% c(issue_sphase, issue_qc, issue_ploidy))]
issues = unique(base::union(base::union(issue_sphase, issue_qc), issue_ploidy))
names = c(pass_names, issues)
stopifnot(length(setdiff(names, colnames(object))) == 0)
ha = rowAnnotation(replicating = ifelse(names %in% issue_sphase, "TRUE", "FALSE"),
                   quality = ifelse(names %in% issue_qc, "TRUE", "FALSE"),
                   ploidy = ifelse(names %in% issue_ploidy, "TRUE", "FALSE"),
                   remove = ifelse(names %in% c(issue_sphase, issue_qc, issue_ploidy), "FALSE", "TRUE"),
                   col = list(replicating = c("FALSE" = "#63ACBE", "TRUE" = "#EE442F"),
                              quality = c("FALSE" = "#63ACBE", "TRUE" = "#EE442F"),
                              ploidy = c("FALSE" = "#63ACBE", "TRUE" = "#EE442F"),
                              remove = c("TRUE" = "#63ACBE", "FALSE" = "#EE442F")),
                   show_legend = c(FALSE, FALSE, FALSE, TRUE))

# optionally, plot raw copy number profiles
#fig_cn = plotCopynumberHeatmap(object[, names], cluster_rows = FALSE, har = ha)
#fig_cn


## select profiles to keep for scUnique analysis
df$keep = !(df$name %in% base::union(base::union(issue_qc, issue_sphase), issue_ploidy))
print("Total kept")
print(prop.table(table(df$keep, df$UID), margin=2))
print(table(df$keep, df$UID))

stopifnot(all(df$rpc[df$keep] > 0 & df$alpha[df$keep] > 0))
stopifnot(all(df$hmm.alpha[df$keep] > 0))
stopifnot(all(!is.na(df$hmm.alpha[df$keep]) & !is.nan(df$hmm.alpha[df$keep])))

ggplot(data = df %>% dplyr::filter(keep)) + geom_histogram(aes(x=ploidy), bins=50) + facet_wrap(~UID+SLX)

## export names of files that pass QC processing to file
print("Writing to file.")
ifelse(!dir.exists(file.path(WORKFLOW_PATH, "results/pass_qc/")), dir.create(file.path(WORKFLOW_PATH, "results/pass_qc/")), FALSE)
readr::write_tsv(df %>% dplyr::filter(keep) %>% dplyr::arrange(name) %>% dplyr::select(UID, SLX, name),
                 path=paste0(WORKFLOW_PATH, "results/pass_qc/pass_", sampleName, ".tsv"), append = FALSE, col_names=FALSE)
