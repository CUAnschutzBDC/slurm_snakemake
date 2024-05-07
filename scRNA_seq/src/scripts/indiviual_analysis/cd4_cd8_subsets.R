library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "combined_celltype"
clusters <- "combined_cluster"
pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6456_PLN"
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

cd4_cells <- c("CD4_Naive", "CD4_TCM", "Treg", "CD4.Proliferating")

cd8_cells <- c("CD8_Naive", "CD8_TCM")

cd4_seurat <- subset(seurat_data, subset = updated_celltype %in% cd4_cells)

cd8_seurat <- subset(seurat_data, subset = updated_celltype %in% cd8_cells)

run_analysis <- function(seurat_object, save_dir){
  # PCA --------------------------------------------------------------------------
  
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- FindVariableFeatures(seurat_object)
  # Update variable features
  # Remove B cell chains
  VariableFeatures(seurat_object) <- VariableFeatures(seurat_object)[!grepl("IG[H|L|K]",
                                                                      VariableFeatures(seurat_object))]
  
  # Remove T cell chains
  VariableFeatures(seurat_object) <- VariableFeatures(seurat_object)[!grepl("TR[A|B][V|J|C]",
                                                                      VariableFeatures(seurat_object))]
  
  seurat_object <- ScaleData(seurat_object)
  
  seurat_object$RNA.weight <- NULL
  seurat_object$RNA_celltype <- NULL
  seurat_object$RNA_celltype_bnd <- NULL
  seurat_object$RNA_celltype_seurat <- NULL
  seurat_object$combined_celltype <- NULL
  seurat_object$combined_celltype_bnd <- NULL
  seurat_object$combined_celltype_seurat <- NULL
  seurat_object$combined_cluster <- NULL
  seurat_object$RNA_cluster <- NULL
  seurat_object$ADT.weight <- NULL
  seurat_object$ADT_celltype <- NULL
  seurat_object$ADT_celltype_bnd <- NULL
  seurat_object$ADT_celltype_seurat <- NULL
  seurat_object$ADT_cluster <- NULL
  
  # Remove previous clustering
  remove_cols <- colnames(seurat_object[[]])[grepl("res\\.[0-9]",
                                                colnames(seurat_object[[]]))]
  
  for (i in remove_cols){
    seurat_object[[i]] <- NULL
  }
  
  
  # Set directories
  ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)
  ifelse(!dir.exists(file.path(save_dir, "images")),
         dir.create(file.path(save_dir, "images")), FALSE)
  ifelse(!dir.exists(file.path(save_dir, "files")),
         dir.create(file.path(save_dir, "files")), FALSE)
  ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
         dir.create(file.path(save_dir, "rda_obj")), FALSE)
  
  
  # PCA of gene expression
  seurat_object <- PCA_dimRed(seurat_object, assay = seurat_assay)
  
  RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                        sample_object = seurat_object)
  
  RNA_pcs <- 20
  
  # UMAP -------------------------------------------------------------------------
  # UMAP of gene expression with final resolution selection
  set.seed(0)
  umap_data <- group_cells(seurat_object, sample, save_dir, nPCs = RNA_pcs,
                           resolution = 0.5, 
                           assay = seurat_assay, HTO = HTO)
  
  seurat_object <- umap_data$object
  
  gene_plots <- umap_data$plots
  
  # Clustifyr --------------------------------------------------------------------
  
  # Information for cell mapping
  all_ref_dir <-
    "/beevol/home/wellskri/Analysis/references/single_cell_references"
  ref_dir <- file.path(all_ref_dir, "pbmc/seurat")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference_l2.csv"),
                      header = TRUE, row.names = 1)
  
  pdf(file.path(save_dir, "images", "celltype_mapping.pdf"))
  
  cluster_res <- name_clusters(seurat_object, ref_mat,
                               save_dir = save_dir,
                               save_name = "celltype_seurat", ADT = FALSE,
                               assay = "RNA",
                               nfeatures = 1000, clusters = "RNA_cluster",
                               plot_type = "rna.umap",
                               features = VariableFeatures(cd4_seurat))
  
  seurat_object <- cluster_res$object
  
  seurat_res_seurat_rna <- cluster_res$RNA
  
  plotDimRed(seurat_object, "RNA_celltype_seurat", plot_type = "rna.umap")
  
  
  dev.off()
  # DE ---------------------------------------------------------------------------
  seurat_object$cluster_celltype <- paste(seurat_object$RNA_cluster,
                                          seurat_object$RNA_celltype_seurat,
                                       sep = "_")
  
  marker_list <- find_write_markers(seurat_object = seurat_object,
                                    meta_col = "cluster_celltype",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "RNA",
                                    save_dir = save_dir)
  
  saveRDS(seurat_object, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
  
  return(seurat_object)
}

cd4_save_dir <- file.path(results_dir, "R_analysis", paste0(sample, "_cd4"))
cd8_save_dir <- file.path(results_dir, "R_analysis", paste0(sample, "_cd8"))


cd4_seurat <- run_analysis(seurat_object = cd4_seurat, save_dir = cd4_save_dir)
cd8_seurat <- run_analysis(seurat_object = cd8_seurat, save_dir = cd8_save_dir)
