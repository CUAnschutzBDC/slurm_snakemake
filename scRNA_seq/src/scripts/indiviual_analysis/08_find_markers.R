library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "combined_celltype"
cell_type2 <- "updated_celltype"
clusters <- "combined_cluster"
pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6553_PLN"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution


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

# Cell type DE -----------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = cell_types,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = cell_types,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# RNA cluster DE ---------------------------------------------------------------

seurat_data$cluster_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[cell_types]][[1]])

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = "cluster_celltype",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = "cluster_celltype",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}


# Cell type DE -----------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = cell_type2,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = cell_type2,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# RNA cluster DE ---------------------------------------------------------------

seurat_data$cluster_updated_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[cell_type2]][[1]])

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = "cluster_updated_celltype",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = "cluster_updated_celltype",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
