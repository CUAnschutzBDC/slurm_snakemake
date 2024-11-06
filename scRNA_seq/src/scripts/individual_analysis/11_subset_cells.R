library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
source(here("src", "scripts", "common_setup.R"))
source(here("src", "scripts", "functions.R"))

subset_name <- "subset_name"
subset_column <- "subset_column"
subset_cell_types <- "subset_cell_types"

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Set directories
save_dir <- file.path(results_dir, "R_analysis", paste(sample, subset_name,
                                                        sep = "_"))

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

# Subset to only endocrine
seurat_data <- subset(seurat_data, subset = get(subset_column) %in%
                        subset_cell_types)

# Clean up
seurat_data <- DietSeurat(seurat_data)
all_celltypes <- colnames(seurat_data[[]])[grepl("RNA_", colnames(seurat_data[[]]))]
all_clusters <- colnames(seurat_data[[]])[grepl("cluster", colnames(seurat_data[[]]))]
all_other <- colnames(seurat_data[[]])[grepl("rna", colnames(seurat_data[[]]))]

if(ADT){
  all_celltypes <- c(
    all_celltypes,
    colnames(seurat_data[[]])[grepl("ADT_", colnames(seurat_data[[]]))],
    colnames(seurat_data[[]])[grepl("combined_", colnames(seurat_data[[]]))]
  )
  
  all_other <- c(
    all_other,
    colnames(seurat_data[[]])[grepl("adt", colnames(seurat_data[[]]))],
    colnames(seurat_data[[]])[grepl("combined", colnames(seurat_data[[]]))]
  )
  
}

for (i in unique(c(all_celltypes, all_clusters, all_other))){
  seurat_data[[i]] <- NULL
}


DefaultAssay(seurat_data) <- "RNA"

seurat_data <- FindVariableFeatures(seurat_data)
seurat_data <- ScaleData(seurat_data)

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                      "seurat_start.rds"))
