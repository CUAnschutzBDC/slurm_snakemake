# Common setup file used in individual and integrated analysis.
#
# This file will set the theme normalization methodm then load in the
# appropriate files and extract fields from the sample_info.csv file
# for use across the pipeline.

library(ggplot2)

# Set theme and select normalization method
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))
normalization_method <- "log" # can be SCT or log

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)

results_dir <- args[[2]]

sample_info <- args[[3]]
sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)
samples_use <- sample_info[sample_info$sample != sample, ]$sample
sample_info <- sample_info[sample_info$sample == sample, ]

# Pull out key fields from sample_info
ADT <- sample_info$ADT
ADT_pca <- sample_info$adt_pca
HTO <- sample_info$HTO
VDJ_T <- sample_info$VDJ_T
VDJ_B <- sample_info$VDJ_B
RNA_pcs <- sample_info$PCs
batch_correction <- sample_info$batch_correction
hash_ident <- sample_info$hash_ident
resolution <- sample_info$resolution
run_adt_umap <- sample_info$adt_umap

# Declare normalization depending on normalization_method
SCT <- normalization_method == "SCT"
seurat_assay <- ifelse(normalization_method == "SCT", "SCT", "RNA")

# Set save directory
save_dir <- file.path(results_dir, "R_analysis", sample)