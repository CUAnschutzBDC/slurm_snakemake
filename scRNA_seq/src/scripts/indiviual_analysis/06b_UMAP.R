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

# UMAP of gene expression with final resolution selection
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = resolution,
                         assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data <- BuildClusterTree(seurat_data, dims = 1:RNA_pcs)
PlotClusterTree(seurat_data)


if(ADT & run_adt_umap){
  if(adt_PCA){
    adt_reduction <- "apca"
  } else{
    adt_reduction <- "pdsb"
    DefaultAssay(seurat_data) <- "ADT"
    ADT_pcs <- nrow(seurat_data)
  }
  
  DefaultAssay(seurat_data) <- "ADT"
  adt_list <- rownames(seurat_data)
  # First set up object
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
                           resolution = 0.2, assay = "ADT", HTO = HTO,
                           reduction = adt_reduction)
  
  seurat_data <- umap_data$object
  
  # Test a variety of resolutions
  seurat_data <- FindClusters(seurat_data, resolution = c(0.1, 0.2, 0.5, 0.8),
                              graph.name = "ADT_snn")
  clustree(seurat_data, prefix = "ADT_snn_res.")
  
  # UMAP of surface protein with selected resolution
  set.seed(0)
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
                           resolution = 0.5, assay = "ADT", HTO = HTO, reduction = adt_reduction)
  
  seurat_data <- umap_data$object
  
  adt_plots <- umap_data$plots
  
  
  plots <- featDistPlot(seurat_data, geneset = adt_list,
                        sep_by = "ADT_cluster", combine = FALSE,
                        assay = "ADT")
  
  # WNN ------------------------------------------------------------------------
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
    seurat_data, reduction.list = list(pca_slot, adt_reduction), 
    dims.list = list(1:RNA_pcs, 1:ADT_pcs),
    modality.weight.name = c(weight_name, "ADT.weight")
  )
  
  
  # Test a wide range of resolutions
  seurat_data <- FindClusters(seurat_data, resolution = c(0.2, 0.5, 0.6, 0.8),
                              graph.name = "wsnn")
  clustree(seurat_data, prefix = "wsnn_res.")
  
  # Run analysis with selected resolution
  seurat_data <- RunUMAP(seurat_data, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_data <- FindClusters(seurat_data, graph.name = "wsnn",
                              algorithm = 3, resolution = 0.6, verbose = FALSE)
  
  seurat_data[["combined_cluster"]] <- Idents(seurat_data)
  col_by_list <- c("combined_cluster", "orig.ident", adt_list)
  if(HTO){
    col_by_list <- c(col_by_list, "HTO_classification")
  }
  save_plot <- file.path(save_dir, "images", "combined_umap.pdf")
  plot_list <- plotDimRed(sample_object = seurat_data,
                          save_plot = save_plot,
                          col_by = col_by_list,
                          plot_type = "wnn.umap")
  VlnPlot(seurat_data, features = "RNA.weight",
          group.by = 'combined_cluster',
          sort = TRUE,
          pt.size = 0.1) +
    NoLegend()
  VlnPlot(seurat_data, features = "ADT.weight",
          group.by = 'combined_cluster',
          sort = TRUE,
          pt.size = 0.1) +
    NoLegend()
  
  DefaultAssay(seurat_data) <- "RNA"
}


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
