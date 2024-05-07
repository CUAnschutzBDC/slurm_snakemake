library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(clustree)
library(harmony)
library(singlecellmethods)
library(batchelor)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))


normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident


ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_adtnorm.rds"))

DefaultAssay(seurat_data) <- "RNA"
# PCA --------------------------------------------------------------------------
VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("IG[H|L|K]",
                                                                      VariableFeatures(seurat_data))]

# PCA of gene expression, weighted by cell number per sample
seurat_data <- singlecellmethods::RunBalancedPCA(obj = seurat_data, 
                                                 weight.by = "orig.ident",
                                                 npcs = 50)

RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                      sample_object = seurat_data,
                      jackstraw = FALSE)

seurat_data$date.processed.for.scSeq <- factor(seurat_data$date.processed.for.scSeq)

all_pc_plots_raster <- plotDimRed(seurat_data, 
                                  col_by = c("percent.mt",
                                             "nFeature_RNA",
                                             "nCount_RNA",
                                             "sample",
                                             "date.processed.for.scSeq"),
                                  plot_type = "pca",
                                  ggrastr = TRUE)

pdf(file.path(save_dir, "images/RNA_pca.pdf"))
print(RNA_plots$pca_loadings)
print(RNA_plots$elbow)
print(all_pc_plots_raster)
dev.off()

if(ADT & run_adt_umap){
  # PCA of surface protein
  seurat_data <- PCA_dimRed(seurat_data, assay = "ADT")
  
  ADT_plots <- plot_PCA(HTO = HTO, assay = "ADT", sample_object = seurat_data)
  
  pdf(file.path(save_dir, "images/ADT_pca.pdf"))
  plot(ADT_plots)
  dev.off()
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
