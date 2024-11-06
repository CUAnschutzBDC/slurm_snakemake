library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)


source(here("src/scripts/common_setup.R"))

# Set to TRUE if you needed to set cell types manually, also
# name the manual cell types "final_celltype"
manual_cell_types <- FALSE

# First, determine the clusters
if (ADT & run_adt_umap) {
  clusters <- "combined_cluster"
} else {
  clusters <- "RNA_cluster"
}

# Then, determine the cell types
if (manual_cell_types) {
  cell_types <- "final_celltype"
} else if (ADT & run_adt_umap) {
  cell_types <- "combined_celltype"
} else {
  cell_types <- "RNA_celltype"
}

pval <- 0.05
logfc <- 0.5

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

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
