library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(clustree)

source(here("src/scripts/common_setup.R"))

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Make base umap
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.8, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

# Test a range of resolutions
seurat_data <- FindClusters(seurat_data, resolution = c(0.2, 0.5, 0.6,
                                                        0.8, 1, 1.2))

clustering_columns <- colnames(seurat_data[[]])[grepl("RNA_snn_res",
                                                      colnames(seurat_data[[]]))]

quality_columns <- c("nCount_RNA", "nFeature_RNA",
                     "percent.mt", "Doublet_finder")

if (ADT){
  quality_columns <- c(quality_columns, "nCount_ADT", "nFeature_ADT")
}

plot_columns <- c(clustering_columns, quality_columns)

all_plots <- plotDimRed(seurat_data, col_by = plot_columns,
                        plot_type = "rna.umap")

# Save these resolutions to a file
pdf(file.path(save_dir, "images", "clustering_resolution.pdf"))
print(clustree(seurat_data))
print(all_plots)

dev.off()

