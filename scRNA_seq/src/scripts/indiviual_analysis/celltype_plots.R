library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(MetBrewer)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

cor_cutoff <- 0.5

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "Npod6456_PLN"

#sample_info <- args[[3]]
sample_info <- here("files/sample_info.tsv")

#results_dir <- args[[2]]
results_dir <- here("results")

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

celltype_colors <- met.brewer(name = "Signac", 
                              n = length(unique(seurat_data$updated_celltype)))

names(celltype_colors) <- unique(seurat_data$updated_celltype)

all_adts <- rownames(GetAssayData(seurat_data, assay = "ADT"))

# Plots to make
# Violin plots of ADTs across new cell types - make sure to name these well
# Some RNA 
# 1. CCR7, TNFRSF4, FOXP3, TCL1A, JCHAIN, BANK1
# nFeature RNA
# nCount RNA
all_violins <- featDistPlot(seurat_data, geneset = all_adts, assay = "ADT",
                            sep_by = "updated_celltype",
                            col_by = "updated_celltype",
                            color = celltype_colors,
                            combine = FALSE)

names(all_violins) <- all_adts

all_violins <- lapply(names(all_violins), function(x){
  save_name <- gsub("ADT-", "", x)
  new_violin <- all_violins[[x]] +
    ggplot2::ggtitle(paste0("Surface protein expression of ", save_name))
  
  return(new_violin)
})

goi <- c("CCR7", "TCL1A", "JCHAIN", "BANK1")


rna_violins <- featDistPlot(seurat_data, geneset = goi, assay = "RNA",
                            sep_by = "updated_celltype",
                            col_by = "updated_celltype",
                            color = celltype_colors,
                            combine = FALSE)

names(rna_violins) <- goi

rna_violins <- lapply(names(rna_violins), function(x){
  new_violin <- rna_violins[[x]] +
    ggplot2::ggtitle(paste0("RNA expression of ", x))
  
  return(new_violin)
})

pdf(file.path(save_dir, "images", "expression_violins_celltype.pdf"))
print(all_violins)
print(rna_violins)
dev.off()

# Violin plots of ADTs across new cell types and clusters - different document
# Some RNA 
# 1. CCR7, TNFRSF4, FOXP3, TCL1A, JCHAIN, BANK1

all_violins <- featDistPlot(seurat_data, geneset = all_adts, assay = "ADT",
                            sep_by = "cluster_updated_celltype",
                            col_by = "updated_celltype",
                            color = celltype_colors,
                            combine = FALSE)

names(all_violins) <- all_adts

all_violins <- lapply(names(all_violins), function(x){
  save_name <- gsub("ADT-", "", x)
  new_violin <- all_violins[[x]] +
    ggplot2::ggtitle(paste0("Surface protein expression of ", save_name))
  
  return(new_violin)
})

goi <- c("CCR7", "TCL1A", "JCHAIN", "BANK1")


rna_violins <- featDistPlot(seurat_data, geneset = goi, assay = "RNA",
                            sep_by = "cluster_updated_celltype",
                            col_by = "updated_celltype",
                            color = celltype_colors,
                            combine = FALSE)

names(rna_violins) <- goi

rna_violins <- lapply(names(rna_violins), function(x){
  new_violin <- rna_violins[[x]] +
    ggplot2::ggtitle(paste0("RNA expression of ", x))
  
  return(new_violin)
})

pdf(file.path(save_dir, "images", "expression_violins_cluster.pdf"))
print(all_violins)
print(rna_violins)
dev.off()

pdf(file.path(save_dir, "images", "umaps.pdf"))
print(plotDimRed(seurat_data, col_by = "updated_celltype",
                 color = celltype_colors, plot_type = "wnn.umap"))

print(plotDimRed(seurat_data, col_by = "cluster_updated_celltype",
                 plot_type = "wnn.umap"))

dev.off()

# 
# plotDimRed(seurat_data, col_by = "ADT-CD8", assay = "ADT",
#            plot_type = "rna.umap")
# 
# plotDimRed(seurat_data, col_by = "ADT-CD8", assay = "ADT",
#            plot_type = "wnn.umap")
# 
# plotDimRed(seurat_data, col_by = "combined_cluster",
#            plot_type = "wnn.umap")
# 
# plotDimRed(seurat_data, col_by = "combined_cluster",
#            plot_type = "wnn.umap", highlight_group = TRUE,
#            group = 14, meta_data_col = "combined_cluster")
# 
# plotDimRed(seurat_data, col_by = "combined_cluster",
#            plot_type = "wnn.umap", highlight_group = TRUE,
#            group = 15, meta_data_col = "combined_cluster")
# 
# 
# plotDimRed(seurat_data, col_by = "percent.mt",
#            plot_type = "wnn.umap")
# 
# 
# plotDimRed(seurat_data, col_by = "nFeature_RNA",
#            plot_type = "wnn.umap")
# 
# plotDimRed(seurat_data, col_by = "nCount_RNA",
#            plot_type = "wnn.umap")
# 
# featDistPlot(seurat_data, geneset = "nFeature_RNA", sep_by= "combined_cluster",
#              col_by= "combined_celltype", combine = FALSE)
# 
# featDistPlot(seurat_data, geneset = "nCount_RNA", sep_by= "combined_cluster",
#              col_by= "combined_celltype", combine = FALSE)
# 
# 
# featDistPlot(seurat_data, geneset = "ITGAX", sep_by= "combined_cluster",
#              col_by= "combined_celltype", combine = FALSE)
# 
# plotDimRed(seurat_data, col_by = "ADT-CD3", assay = "ADT",
#            plot_type = "wnn.umap")
# 
# plotDimRed(seurat_data, col_by = "ADT-CD4", assay = "ADT",
#            plot_type = "wnn.umap")
# 
# 
# plotDimRed(seurat_data, col_by = "XIST", assay = "RNA",
#            plot_type = "wnn.umap")
# 
# plotDimRed(seurat_data, col_by = "CD3D", assay = "RNA",
#            plot_type = "wnn.umap")
# 
# plotDimRed(seurat_data, col_by = "ITGAX", assay = "RNA",
#            plot_type = "wnn.umap")
# 
# 
# plotDimRed(seurat_data, col_by = "CD4", assay = "RNA",
#            plot_type = "wnn.umap")
# 
# 
# plotDimRed(seurat_data, col_by = "CD8A", assay = "RNA",
#            plot_type = "wnn.umap")
# 
