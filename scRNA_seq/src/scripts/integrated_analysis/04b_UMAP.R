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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# UMAP -------------------------------------------------------------------------
seurat_data$RNA_cluster <- NULL

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}

## Cluster --------------------------------------------------------------------

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = resolution, assay = seurat_assay,
                         HTO = HTO,
                         reduction = batch_correction)
seurat_data <- umap_data$object

plot_type <- paste0(batch_correction, ".umap")

seurat_data$corrected_cluster <- seurat_data$RNA_cluster

if(ADT & run_adt_umap){
  # UMAP of surface protein
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
                           resolution = 0.6, assay = "ADT", HTO = TRUE)
  
  seurat_data <- umap_data$object
  
  adt_plots <- umap_data$plots
  
  
  # UMAP of combined
  # Identify multimodal neighbors. These will be stored in the neighbors slot, 
  # and can be accessed using bm[['weighted.nn']]
  # The WNN graph can be accessed at bm[["wknn"]], 
  # and the SNN graph used for clustering at bm[["wsnn"]]
  # Cell-specific modality weights can be accessed at bm$RNA.weight
  if(SCT){
    pca_slot <- "sctpca"
    weight_name <- "SCT.weight"
  } else{
    pca_slot <- "pca"
    weight_name <- "RNA.weight"
  }
  seurat_data <- FindMultiModalNeighbors(
    seurat_data, reduction.list = list(pca_slot, "apca"), 
    dims.list = list(1:RNA_pcs, 1:ADT_pcs),
    modality.weight.name = c(weight_name, "ADT.weight")
  )
  
  
  seurat_data <- RunUMAP(seurat_data, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_data <- FindClusters(seurat_data, graph.name = "wsnn",
                              algorithm = 3, resolution = 2, verbose = FALSE)
  
  seurat_data[["combined_cluster"]] <- Idents(seurat_data)
  col_by_list <- c("combined_cluster", "orig.ident")
  if(HTO){
    col_by_list <- c(col_by_list, "HTO_classification")
  }
  save_plot <- file.path(save_dir, "images/combined_umap.pdf")
  plot_list <- plotDimRed(sample_object = seurat_data,
                          save_plot = NULL,
                          col_by = col_by_list, return_plot = TRUE,
                          plot_type = "wnn.umap")
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))
