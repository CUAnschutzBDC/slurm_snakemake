library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source(here("src/scripts/functions.R"))

remove_ambience <- FALSE

normalization_method <- "log" # can be SCT or log

cor_cutoff <- 0.5

all_ref_dir <-
  "/pl/active/Anschutz_BDC/resources/single_cell_references"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

# ADT_pca <- sample_info$adt_pca
# run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution
batch_correction <- sample_info$batch_correction

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Set directories
save_dir <- file.path(results_dir, "R_analysis", paste0(sample, "_EP"))

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

# Subset to only endocrine
seurat_data <- subset(seurat_data, subset = RNA_celltype %in%
                        c("EP1"))

# Clean up
seurat_data <- DietSeurat(seurat_data)
all_celltypes <- colnames(seurat_data[[]])[grepl("RNA_", colnames(seurat_data[[]]))]
all_clusters <- colnames(seurat_data[[]])[grepl("cluster", colnames(seurat_data[[]]))]
all_other <- colnames(seurat_data[[]])[grepl("rna", colnames(seurat_data[[]]))]

for (i in unique(c(all_celltypes, all_clusters, all_other))){
  seurat_data[[i]] <- NULL
}


seurat_data <- FindVariableFeatures(seurat_data)
seurat_data <- ScaleData(seurat_data)

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                      "seurat_norm.rds"))
