# Add your samples in here if you think the cell types need to be updated,
# if your sample isn't in the list, udpates will be skipped and the 
# "updated_celltype" will just be set to the value in the "cell_types"

library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
source(here("src", "scripts", "common_setup.R"))

if(ADT){
  cell_types <- "combined_celltype"
  clusters <- "combined_cluster"
} else {
  cell_types <- "RNA_celltype"
  clusters <- "RNA_cluster"
}

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
pdf(file.path(save_dir, "images", "cell_type_plots.pdf"))

plot(featDistPlot(seurat_data, 
             geneset = c("JCHAIN", "XBP1", "IRF4", "nFeature_RNA"),
             sep_by= clusters,
             col_by= clusters, combine = TRUE))

# From https://www.science.org/doi/full/10.1126/science.abf1970
plot_list <- c("CD4", "CD40LG", # CD4+ cells
               "CD8A", "CD8B", # CD8+ cells
               "CCR7", # Naive
               "IL7R", "TNFRSF4", # EM
               "RTKN2", "FOXP3", # CD4 reg
               "PRF1", 
               "GZMH", # CD8 GZMH
               "GZMB", 
               "GZMK", # CD8 GZMK
               "KLRB1", # CD8 MAIT
               "GNLY", "NKG7", "FCGR3A", #NK
               "TCL1A", # B naive
               "BANK1", # B mem
               "MZB1", # B plasma
               "FCRL5",
               "CD19", "MS4A1", "CD79A", # B other
               "SOX4")

print(featDistPlot(seurat_data, geneset = "JCHAIN", sep_by= "combined_cluster",
             col_by= "combined_celltype", combine = FALSE))

all_plots <- featDistPlot(seurat_data, geneset = plot_list, sep_by= "combined_cluster",
                          col_by= "combined_celltype", combine = FALSE)

print(all_plots)

if (ADT){
  DefaultAssay(seurat_data) <- "SCAR_ADT"
  all_adts <- rownames(seurat_data)
  DefaultAssay(seurat_data) <- seurat_assay
  
  adt_plots <- featDistPlot(seurat_data, geneset = all_adts, sep_by= "combined_cluster",
                            col_by= "combined_celltype", combine = FALSE,
                            assay = "SCAR_ADT")  
}


# CD4 - 16, 5, 4, 9, 0, 2, 15 (low)
# CD40LG - 4, 5 (low), 3 (low), 9, 0, 2
# CD8A - 1, 12, 6
# CD8B - 1, 12, 6
# CCR7 - 0, 1, 10 (low), 12, 15, 2, 3, 4 (low), 5 (low), 9
# IL7R - 0, 1, 11, 12, 13, 14, 15, 16, 2, 3, 4, 5, 6, 9
# TNFRSF4 - 5, 4
# RTKN2 - 5, 14
# FOXP3 - 16, 5
# PRF1 - 13
# GZMH 16
# GZMB 13
# GZMK - 13, 6, 16, 11
# GNLY - 13
# NKG7 - 13, 6
# FCGR3A - 13, 6
# TCL1A - 7
# BANK1 - 7, 10, 8
# MZB1 - 14
# FCRL5 - 14, 10 (low)
# CD19 - 15, 7, 10, 8
# MS4A1 - 15, 7, 10, 8
# CD79A - 7, 10, 8
# SOX4 - 13


cluster_colors <- colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9,name = "Set1")
  )(length(unique(seurat_data[[clusters]][[1]])))

names(cluster_colors) <- unique(seurat_data[[clusters]][[1]])

celltype_colors <- colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9,name = "Set1")
)(length(unique(seurat_data[[cell_types]][[1]])))

names(celltype_colors) <- unique(seurat_data[[cell_types]][[1]])

coloring <- list(clusters = cluster_colors,
                 cell_types = celltype_colors)


meta_data <- seurat_data[[]] %>%
  dplyr::select(!!sym(clusters), !!sym(cell_types)) %>%
  dplyr::mutate(cluster_celltype = paste(!!sym(clusters), !!sym(cell_types),
                                         sep = "_")) %>%
  dplyr::distinct()

rownames(meta_data) <- meta_data$cluster_celltype

seurat_data$cluster_celltype <- paste(seurat_data[[clusters]][[1]],
                                      seurat_data[[cell_types]][[1]],
                                      sep = "_")

plot_heatmap(seurat_object = seurat_data, gene_list = plot_list,
             meta_df = meta_data,
             meta_col = "cluster_celltype", cluster_cols = TRUE,
             average_expression = TRUE,
             plot_meta_col = FALSE, color_list = coloring)

dev.off()

if(sample == "Npod6456_PLN"){
  # Final cell type mapping
  # CLuster 0 = CD4_Naive - based on clustifyr and ADTs and RNA (CCR7)
  # Cluster 1 = CD8 Niave - based on clustifyr, ADTs, and RNA (CCR7)
  # Cluster 2 = CD4_Naive - based on clustifyr and ADTs and RNA (CCR7)
  # Cluster 3 = Low quality - based on number of features, DE genes, number of counts
  # Cluster 4 = CD4_TCM - based on clustifyr and ADTs and RNA (TNFRSF4)
  # Cluster 5 = Treg - based on clustifyr and ADTs and RNA (FOXP3)
  # Cluster 6 = CD8_TCM - based on clustifyr and ADTs
  # Cluster 7 = Naive 1 - based on clustifyr and ADTs and RNA (TCL1A)
  # Cluster 8 = B memory - based on clustifyr, ADTs, and RNA (BANK1)
  # Cluster 9 = CD4_Naive - based on clustifyr and ADTs and RNA (CCR7)
  # Cluster 10 = B memory - based on clustifyr, ADTs, and RNA (BANK1)
  # Cluster 11 = Unknown
  # Cluster 12 = CD8_Naive - based on clustifyr and ADTs and RNA (CCR7)
  # Cluster 13 = NK
  # Cluster 14 = plasmablast - based on gene expression (JCHAIN) and DE
  # Cluster 15 = doublet
  updated_celltype = c("0" = "CD4_Naive",
                       "1" = "CD8_Naive",
                       "2" = "CD4_Naive",
                       "3" = "low_quality",
                       "4" = "CD4_TCM",
                       "5" = "Treg",
                       "6" = "CD8_TCM",
                       "7" = "B_Naive",
                       "8" = "B_Memory",
                       "9" = "CD4_Naive",
                       "10" = "B_Memory",
                       "11" = "Unknown",
                       "12" = "CD8_Naive",
                       "13" = "NK",
                       "14" = "Plasmablast",
                       "15" = "Doublet",
                       "16" = "CD4_Unknown")
} else if (sample == "Npod6553_PLN"){
  # Final cell type mapping
  # Cluster 0 - plasmablast (JCHAIN, clustifyr)
  # Cluster 1 - Treg (CD3, CD4)
  # Cluster 2 - Treg (CD3, CD4)
  # Cluster 3 - T cell cd8 naive (CD3, CD8, CCR7, clustifyr)
  # Cluster 4 - doublet
  # Cluster 5 - B cell memory (CXCR5, CD21, IGM, IGD, BANK1)
  # Cluster 6 - other (low CD3, CD4, CD8)
  # Cluster 7 - T cell cd8 (CD3, CD8)
  # Cluster 8 - other
  
  updated_celltype = c("0" = "plasmablast",
                       "1" = "Treg",
                       "2" = "Treg",
                       "3" = "CD8_Naive",
                       "4" = "Doublet",
                       "5" = "B_memory",
                       "6" = "Other",
                       "7" = "CD8_TCM",
                       "8" = "ILC")
} else {
  updated_celltype <- seurat_data[[cell_types]][[1]]
  names(updated_celltype) <- seurat_data[[clusers]][[1]]
  updated_celltype <- updated_celltype[!duplicated(updated_celltype)]
}



seurat_data$updated_celltype <- updated_celltype[seurat_data[[clusters]][[1]]]


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
